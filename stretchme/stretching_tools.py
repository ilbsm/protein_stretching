""" The miscellanea functions for stretching package"""

import pandas as pd
import numpy as np
from .default_parameters import default_parameters
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit, minimize, fmin
from scipy.special import erf
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
from colorsys import hsv_to_rgb
from scipy.stats import cauchy, ks_2samp
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from lmfit import Model


# preprocessing
def pack_parameters(filename, sheet_name=0, linker=None, unit='nm', speed=1, residues_distance=0.365, source='theory',
                    low_force_cutoff=0.1, plot_columns=4, initial_guess=None, separator=None):

    """ Filtering and packing the parameters for a clearer view."""
    # TODO add filters of the input parameters
    parameters = {
        'filename': filename,
        'sheet_name': sheet_name,
        'linker': linker,
        'unit': unit,
        'speed': speed,
        'source': source,
        'residues_distance': residues_distance,
        'plot_columns': plot_columns,
        'separator': separator,
        'low_force_cutoff': low_force_cutoff
    }
    if initial_guess:
        parameters['initial_guess'] = initial_guess
    else:
        parameters['initial_guess'] = default_parameters['initial_guess'][parameters['source']]
    return parameters


def running_average(x, y, window=None):
    if not window:
        window = max(int(len(x) / 100), 8)
    x_smooth = np.convolve(x, np.ones((window,)) / window, mode='valid')
    y_smooth = np.convolve(y, np.ones((window,)) / window, mode='valid')
    return x_smooth, y_smooth


# calculations
def gauss(x, height, mean, width):
    return height * np.exp(-((x - mean) ** 2) / (2 * width))


def integrate_gauss(force, mean, width):
    return 0.5 * (1 - erf((force-mean)/(np.sqrt(width * 2))))


def invert_wlc_np(forces, p, k=0):
    return np.array([invert_wlc(f, p, k) for f in forces])


def invert_wlc(force, p, k=0):
    if k == 0:
        coefs = [1, -(2.25 + force / p), (1.5 + 2 * force / p), -force / p]
    else:
        coefs = [1,
                 -(2.25 + force * (3*k + 1 / p)),
                 (3 * (k ** 2) + 2 * (k / p)) * force ** 2 + ((4.5 * k) + (2 / p)) * force + 1.5,
                 -force * (((k**3) + ((k**2) / p)) * (force ** 2) +
                           (2.25 * (k**2) + 2 * (k / p)) * force +
                           ((1.5 * k) + (1 / p)))]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    if k == 0:
        result = result[result < 1]
    return min(result)


def ewlc(d, length, p, k):
    if k == 0:
        if d < 0.99 * length:
            return p * (0.25 / ((1 - d / length) ** 2) - 0.25 + d / length)
        else:
            return 999
    else:
        x = d/length
        return reduced_ewlc(x, p, k)


def reduced_ewlc(x, p, k):
    coefs = [-(k**3) - (k**2)/p,
             -(2.25 * (k**2)) - 2 * k/p + x * (3 * (k**2) + 2 * (k/p)),
             -(1.5 * k) - 1/p + x * (4.5 * k + 2/p) - x**2 * (3 * k + 1/p),
             1.5 * x - 2.25 * x**2 + x**3]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    return max(result)


def reduced_ewlc_np(x, p, k):
    return np.array([reduced_ewlc(point, p, k) for point in x])


def ewlc_np(d, p, k, l):
    return np.array([reduced_ewlc(point/l, p, k) for point in d])


def wlc_np(d, p, l):
    return np.array([reduced_wlc(point/l, p) for point in d])


def reduced_wlc(x, p):
    return p * (0.25 / ((1 - x) ** 2) - 0.25 + x)


def get_d_dna(p_dna, l_dna, k_dna, f_space):
    if p_dna > 0:
        return l_dna * np.array([invert_wlc(f, p_dna, k_dna) for f in f_space])
    else:
        return np.zeros(len(f_space))


def get_force_load(force, speed, l_prot, p_prot, k_prot=None):
    numerator = 2 * (l_prot/p_prot) * (1 + force/p_prot)
    denominator = 3 + 5 * force/p_prot + 8 * (force/p_prot)**(5/2)
    factor = numerator/denominator
    if k_prot:
        factor += 1/k_prot
    return speed/factor


def dhs_feat_05(force, x, t0, g):
    return t0 / (1 - 0.5 * x/g * force) * np.exp(-g*(1-(1-0.5 * x/g * force)**2))


def dhs_feat_066(force, x, t0, g):
    return t0 / (1 - 2 * x/g * force/3)**(1/2) * np.exp(-g*(1-(1-2 * x/g * force/3)**(3/2)))


def dhs_feat_1(force, x, t0):
    return t0 * np.exp(-x*force)


def dhs_feat(f_space, lifetime):
    coefficients = {}
    # # v = 1
    # p0 = (1, lifetime[0])
    # try:
    #     popt, pcov = curve_fit(dhs_feat_1, f_space, lifetime, p0=p0)
    #     coefficients['1'] = {'x': popt[0], 't0': popt[1], 'covariance': pcov}
    # except RuntimeError:
    #     coefficients['1'] = None

    # # v = 2/3
    # p0 = (0.1, lifetime[0], 10)
    # try:
    #     popt, pcov = curve_fit(dhs_feat_066, f_space, lifetime, p0=p0)
    #     coefficients['2/3'] = {'x': popt[0], 't0': popt[1], 'g': popt[2], 'covariance': pcov}
    # except RuntimeError:
    #     coefficients['2/3'] = None

    # v = 1/2
    p0 = (1, lifetime[0], 1)
    try:
        popt, pcov = curve_fit(dhs_feat_05, f_space, lifetime, p0=p0)
        coefficients['1/2'] = {'x': popt[0], 't0': popt[1], 'g': popt[2], 'covariance': pcov}
    except RuntimeError:
        coefficients['1/2'] = None

    return coefficients


# plotting and colors
colors = [mcolors.CSS4_COLORS['red'],
          mcolors.CSS4_COLORS['green'],
          mcolors.CSS4_COLORS['blue'],
          mcolors.CSS4_COLORS['yellow'],
          mcolors.CSS4_COLORS['cyan'],
          mcolors.CSS4_COLORS['orange'],
          mcolors.CSS4_COLORS['purple'],
          mcolors.CSS4_COLORS['lime'],
          mcolors.CSS4_COLORS['magenta']]


def get_gray_shade(k, max_value):
    values = np.linspace(0, 0.75, max_value)
    return hsv_to_rgb(0, 0, values[k])


def get_color(k, max_value):
    if k >= max_value:
        raise ValueError("The number of color requested exceeded the total number of colors.")
    if max_value > 9:
        values = np.linspace(0, 0.66, max_value)
        return hsv_to_rgb(values[k], 1, 1)
    else:
        return colors[k]


def plot_decomposed_histogram(position, data, l_space, residue_distance):
    k = 0
    # TODO add Cauchy distribution
    for index, row in data[['means', 'widths', 'heights', 'gammas']].iterrows():
        mean, width, height, gamma = tuple(row.to_numpy())
        y_plot = gauss(l_space, height, mean, width)
        y_plot_cauchy = cauchy.pdf(l_space, mean, gamma)
        residues = 1 + int(mean / residue_distance)
        label = "L= " + str(round(mean, 3)) + ' (' + str(residues) + ' AA)'
        position.plot(l_space, y_plot, linestyle='--', label=label, color=get_color(k, len(data)))
        k += 1
    return


def plot_trace_fits(position, coefficients, f_space, residue_distance):
    d_dna = get_d_dna(coefficients.get('p_dna', 0), coefficients.get('l_dna', 0),
                      coefficients.get('k_dna', None), f_space)
    x_prot = np.array([invert_wlc(f, coefficients.get('p_prot', 0), coefficients.get('k_prot', None))
                       for f in f_space])
    k = 0
    for index, row in coefficients['l_prot'].iterrows():
        l_prot = row['means']
        residues = 1 + int(l_prot / residue_distance)
        d_prot = l_prot * x_prot
        d_plot = d_dna + d_prot
        label = 'Fit L=' + str(round(l_prot, 3)) + ' (' + str(residues) + ' AA)'
        position.plot(d_plot, f_space, label=label, color=get_color(k, len(coefficients['l_prot'])))
        k += 1
    return


# postprocessing
def decompose_histogram(data, significance, states=None):
    # TODO take care of example with 1 observation...
    if len(data) == 1:
        parameters = pd.DataFrame({'means': np.array([float(data)]), 'widths': np.array([0]), 'heights': np.array(1)})
        boundaries = [float(data), float(data)]
        return parameters, boundaries
    # finding the number of components
    x = np.expand_dims(data, 1)
    kde = KernelDensity().fit(x)
    estimator = np.linspace(min(data), max(data), 1001)
    kde_est = np.exp(kde.score_samples(estimator.reshape(-1, 1)))
    maximas = [estimator[_] for _ in argrelextrema(kde_est, np.greater)[0] if kde_est[_] > significance]
    if not states:
        states = len(maximas)

    # finding the range
    # TODO clean it up
    # TODO take care for the case with len(maximas) == 0
    min_list = [estimator[_] for _ in range(len(estimator)) if kde_est[_] < significance and estimator[_] < maximas[0]]
    max_list = [estimator[_] for _ in range(len(estimator)) if kde_est[_] < significance and estimator[_] > maximas[-1]]
    if len(min_list) > 0:
        min_boundary = max(min_list)
    else:
        min_boundary = min(estimator)
    if len(max_list) > 0:
        max_boundary = min(max_list)
    else:
        max_boundary = max(estimator)

    # fitting Gaussians
    x = np.expand_dims(data[(data > min_boundary) & (data < max_boundary)], 1)
    gmm = GaussianMixture(n_components=states)  # gmm for two components
    gmm.fit(x)
    # finding parameters
    ranges = [[0, 35], [35, 60], [60, 120]]
    cauchy_mus = []
    cauchy_gammas = []
    pvalues = []
    for x in gmm.means_:
        for r in ranges:
            if r[0] < x[0] < r[1]:
                data_hist = data[(data >= r[0]) & (data <= r[1])]
                mu, gamma = cauchy.fit(data_hist)
                cauchy_gammas.append(gamma)
                cauchy_mus.append(mu)
                ensamble = cauchy.rvs(size=len(data_hist), loc=mu, scale=gamma)
                pvalues.append(ks_2samp(ensamble, data_hist).pvalue)
    parameters = pd.DataFrame({'means': np.array([x[0] for x in gmm.means_]),
                               'widths': np.array([x[0][0] for x in gmm.covariances_]),
                               'heights': np.array([np.exp(gmm.score_samples(np.array(u).reshape(-1, 1)))[0]
                                                    for u in [x[0] for x in gmm.means_]]),
                               'gammas': np.array(cauchy_gammas),
                               'cauchy_mus': np.array(cauchy_mus),
                               'pvalues': np.array(pvalues)})
    parameters = parameters.sort_values(by=['means'])
    boundaries = [min_boundary, max_boundary]
    # print(parameters)
    # print(parameters['gammas'].min(), parameters['pvalues'].max())
    return parameters, boundaries


# def implicite_elastic_wlc(data, l, k, p):
#     d = np.array(data['d'])
#     f = np.array(data['F'])
#     x = d / l
#     coefficients = [1,
#                     -f * (3 / k + 1 / p) - 9 / 4,
#                     f ** 2 * (3 / (k ** 2) + 2 / (k * p)) + f * (9 / (2 * k) + 2 / p) + 3 / 2,
#                     -f ** 3 * (1 / (k ** 3) + 1 / (p * k ** 2)) - f ** 2 * (9 / (4 * k ** 2) + 2 / (k * p)) - f * (
#                                 3 / (2 * k) + 1 / p)]
#     return x ** 3 * coefficients[0] + x ** 2 * coefficients[1] + x * coefficients[2] + coefficients[3]
#
#
# def implicite_elastic_wlc_amplitude(data, l, k, p):
#     return np.abs(implicite_elastic_wlc(data, l, k, p))
#
#
# def minimize_kp(df, length, init_k, init_p):
#     ydata = np.zeros(len(df))
#     popt, pcov = curve_fit(implicite_elastic_wlc_amplitude, df, ydata, bounds=(0, np.inf), p0=(length, init_k, init_p))
#     return popt, pcov

def find_state_boundaries(smooth_data, parameters, p, k, max_distance=0.3):
    begs = []
    ends = []
    last_end = 0
    for index, row in parameters.iterrows():
        l_prot = row['means']
        smooth_data['state_' + str(index)] = np.array([ewlc(d, l_prot, p, k) for d in list(smooth_data['d'])])
        data_close = smooth_data[abs(smooth_data['F'] - smooth_data['state_' + str(index)]) <= max_distance]
        data_close = data_close[data_close['d'] > last_end]['d']
        begs.append(data_close.min())
        ends.append(data_close.max())
        if ends[-1] != np.NaN:
            last_end = ends[-1]
    return begs, ends, smooth_data


def minimize_pk(data, data_smoothed, p, k, significance=0.01, max_distance=0.3):
    p_new = [p, k, 0]
    known_results = []
    lmfitmodel = Model(ewlc_np)
    params = lmfitmodel.make_params()
    params['p'].min = 0.1
    params['k'].min = 0.0
    params['l'].min = 0.1
    while True:
        # finding the approximate contour length
        known_results.append(p_new)
        hist_values = data['d']/np.array([invert_wlc(f, p_new[0], p_new[1]) for f in data['F']])
        parameters, boundaries = decompose_histogram(np.array(hist_values), significance)
        means = [round(_, 2) for _ in parameters['means'].values]
        if means in known_results:
            break
        else:
            known_results.append(means)

        # finding ranges of states
        smoothed = pd.DataFrame({'d': data_smoothed['d'], 'F': data_smoothed['F']})
        begs, ends, smoothed = find_state_boundaries(smoothed, parameters, p, k, max_distance)
        parameters['begs'] = begs
        parameters['ends'] = ends

        # selecting the last state
        row = parameters.dropna().sort_values(by='means', ascending=False).head(1)
        beg, end, l_prot = row[['begs', 'ends', 'means']].to_numpy()[0]
        fit_data = pd.DataFrame({'d': data_smoothed[data_smoothed['d'].between(beg, end)]['d'],
                                 'F': data_smoothed[data_smoothed['d'].between(beg, end)]['F']})
        fit_data = fit_data.dropna()

        # fitting
        if p_new[2] >= 0.1:
            p_old = p_new
        else:
            p_old = [p_new[0], p_new[1], float(row['means'])]
        result = lmfitmodel.fit(fit_data['F'], d=fit_data['d'], p=p_old[0], k=p_old[1], l=p_old[2])
        # print(result.fit_report())
        p_new = [round(result.params[key].value, 4) for key in ['p', 'k', 'l']]
    return p_new[:-1]


def implicit_elastic_wlc(data, p, k):
    result = pd.DataFrame(columns=['0', '1', '2', '3', 'sum'])
    result['3'] = data['x'] ** 3
    result['2'] = (-1) * data['x'] ** 2 * (2.25 + data['F'] * (3 / k + 1 / p))
    result['1'] = data['x'] * (data['F'] ** 2 * (3 / (k**2) + 2 / (k * p)) + data['F'] * (4.5 / k + 2 / p) + 1.5)
    result['0'] = data['F'] * ((1 / (k**3) + 1 / (p * k**2)) * data['F'] ** 2 +
                               (2.25 / (k**2) + 2 / (k * p)) * data['F'] + (1.5 / k + 1 / p))
    result['sum'] = np.abs(result['3'] + result['2'] + result['1'] + result['0'])
    return result['sum']


def guess_states_number(hist_values, significance=0.01):
    # TODO take care for the case with len(maximas) == 0

    # finding maximas
    x = np.expand_dims(hist_values, 1)
    kde = KernelDensity().fit(x)
    estimator = np.linspace(min(hist_values), max(hist_values), 1001)
    kde_est = np.exp(kde.score_samples(estimator.reshape(-1, 1)))
    maximas = [estimator[_] for _ in argrelextrema(kde_est, np.greater)[0] if kde_est[_] > significance]

    # finding support
    significant = [estimator[_] for _ in range(len(estimator)) if kde_est[_] >= significance] + maximas
    if significant:
        support = [min(significant), max(significant)]
    else:
        support = [min(estimator), max(estimator)]
    return maximas, support


def decompose_histogram_cauchy(hist_values, states=None, significance=0.01):
    if not states:
        maximas, support = guess_states_number(hist_values, significance=significance)
        states = max(len(maximas), 1)
    else:
        # TODO check if this is enough
        significant = [x for x in hist_values if x >= significance]
        support = [min(significant), max(significant)]

    # fitting Gaussians
    trimmed_data = hist_values[(hist_values > support[0]) & (hist_values < support[1])]
    x = np.expand_dims(trimmed_data, 1)
    gmm = GaussianMixture(n_components=states)
    gmm.fit(x)

    # defining parameters
    parameters = pd.DataFrame({'means': np.array([x[0] for x in gmm.means_]),
                               'widths': np.array([x[0][0] for x in gmm.covariances_])})
    parameters['heights'] = 1/(parameters['widths'] * np.sqrt(2*np.pi))
    parameters = parameters.sort_values(by=['means'])

    states = pd.DataFrame({'d': trimmed_data})
    states_names = []

    for ind, row in parameters.iterrows():
        states_names.append('state_' + str(ind))
        mean, width, height = row[['means', 'widths', 'heights']].values
        states[states_names[-1]] = gauss(trimmed_data, height, mean, width)
    parameters['states_names'] = np.array(states_names)

    # assigning the best matching state
    states['state'] = states[states_names].idxmax(axis=1)

    # fitting Cauchy
    cauchy_means = []
    cauchy_gammas = []
    pvalues = []
    boundary_means = [0] + list(parameters['means'].values) + [max(trimmed_data)]
    for k in range(len(states_names)):
        state = states_names[k]
        bounds = [boundary_means[k], boundary_means[k+2]]
        matching = states[(states['state'] == state) & (states['d'].between(bounds[0], bounds[1]))]['d']
        mu, gamma = cauchy.fit(matching.to_numpy())
        ensamble = cauchy.rvs(size=len(matching), loc=mu, scale=gamma)
        cauchy_means.append(mu)
        cauchy_gammas.append(gamma)
        pvalues.append(ks_2samp(ensamble, matching).pvalue)
    parameters['cauchy_means'] = np.array(cauchy_means)
    parameters['cauchy_gammas'] = np.array(cauchy_gammas)
    parameters['pvalues'] = np.array(pvalues)
    return parameters


def fit_error(x, data, states=None, logger=None):
    if logger:
        logger.info(str(x))
    x_prot = invert_wlc_np(data['F'], x[0], x[1])
    hist_values = data['d']/x_prot
    parameters = decompose_histogram_cauchy(hist_values, states=states)
    return parameters['cauchy_gammas'].min() - parameters['pvalues'].max()


def fit_pk_linker_none(data, p, k, states=None, logger=None):
    if logger:
        logger.info("Fitting the parameters\np\tk")
    x_opt = minimize(fit_error, x0=np.array([p, k]), args=(data, states, logger), bounds=((0.1, 100), (0, 1)))
    p_opt, k_opt = x_opt.x
    return p_opt, k_opt
