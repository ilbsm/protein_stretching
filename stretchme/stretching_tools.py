""" The miscellanea functions for stretching package"""

import pandas as pd
import numpy as np
import logging
import matplotlib.colors as mcolors
from os import path
from io import StringIO
from scipy.optimize import curve_fit, minimize
from scipy.special import erf
from scipy.stats import cauchy, ks_2samp, skew
from scipy.signal import argrelextrema
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from colorsys import hsv_to_rgb
from .default_parameters import default_parameters
from copy import deepcopy
from matplotlib import pyplot as plt


# preprocessing
def pack_parameters(filename, sheet_name=0, linker=None, unit='nm', speed=1, residues_distance=0.365, source='theory',
                    low_force_cutoff=0.1, plot_columns=4, initial_guess=None, method='marko-siggia', separator=None):

    """ Filtering and packing the parameters for a clearer view."""
    # TODO add filters of the input parameters, include all the paramters
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
        'low_force_cutoff': low_force_cutoff,
        'method': method
    }
    if initial_guess:
        parameters['initial_guess'] = initial_guess
    else:
        parameters['initial_guess'] = deepcopy(default_parameters['initial_guess'][parameters['source']])
    return parameters


def running_average(x, y, window=None):
    if not window:
        window = max(int(len(x) / 100), 8)
    x_smooth = np.convolve(x, np.ones((window,)) / window, mode='valid')
    y_smooth = np.convolve(y, np.ones((window,)) / window, mode='valid')
    return x_smooth, y_smooth


# reading
def read_dataframe(input_data, cases=None, columns=None):
    if columns:
        return input_data[columns]
    elif cases:
        allowed = [str(_) for _ in cases]
        colnames = [name for name in list(input_data) if name.strip('dF_') in allowed]
        return input_data[colnames]
    else:
        return input_data


def read_from_file(input_data, cases, columns, parameters, name, debug):
    if not input_data:
        return pd.DataFrame(columns=['d', 'F'])
    elif isinstance(input_data, str) and path.isfile(input_data) and '.xls' in input_data:
        if debug:
            logger = set_logger(name)
            logger.info("Reading the file " + input_data + " as a xls file.")
            close_logs(logger)
        return read_excel(input_data, cases, columns, parameters)
    elif isinstance(input_data, str) and path.isfile(input_data):
        if debug:
            logger = set_logger(name)
            logger.info("Reading the file " + input_data + " as a csv file.")
            close_logs(logger)
        return read_csv(input_data, cases, columns, parameters)
    else:
        raise NotImplementedError("Either file not found, or data not given as Pandas Dataframe.")


def read_excel(input_data, cases, columns, parameters):
    data = pd.read_excel(input_data, sheet_name=parameters['sheet_name'])
    return read_dataframe(data, cases=cases, columns=columns)


def read_csv(input_data, cases, columns, parameters):
    separator = parameters['separator']
    with open(input_data, 'r') as myfile:
        content = myfile.read().split("#")[-1].strip()
    if separator != ' ':
        data = pd.read_csv(StringIO(content), sep=separator, escapechar='#')
    else:
        data = pd.read_csv(StringIO(content), delim_whitespace=True, escapechar='#')
    return read_dataframe(data, cases=cases, columns=columns)


# logging
def set_logger(name, mode='a'):
    logging.basicConfig(filename=name + '.log', filemode=mode, level=logging.DEBUG,
                        format='%(asctime)s %(message)s')
    return logging.getLogger()


def close_logs(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)
    return


# calculations
def gauss(x, height, mean, width):
    return height * np.exp(-((x - mean) ** 2) / (2 * width))


def integrate_gauss(force, mean, width):
    return 0.5 * (1 - erf((force-mean)/(np.sqrt(width * 2))))


# calculating F(d)
def marko_siggia(d, length, p, k=0):
    if k == 0 and d > 0.99 * length:
        return 999
    elif k == 0:
        return p * (0.25 / ((1 - d / length) ** 2) - 0.25 + d / length)
    else:
        x = d/length
        coefs = [-(k ** 3) - (k ** 2) / p,
                 -(2.25 * (k ** 2)) - 2 * k / p + x * (3 * (k ** 2) + 2 * (k / p)),
                 -(1.5 * k) - 1 / p + x * (4.5 * k + 2 / p) - x ** 2 * (3 * k + 1 / p),
                 1.5 * x - 2.25 * x ** 2 + x ** 3]
        result = np.roots(coefs)
        result = np.real(result[np.isreal(result)])
        result = result[result > 0]
        return max(result)


def stretch_adjusted_wlc(d, length, p, k=0, residues_distance=None):
    if not residues_distance:
        residues_distance = default_parameters['residues_distance']
    x = d / length
    g = k * residues_distance
    coefs = [((4 * g**2 * x**2)/p + (4 * g**3 * x**3)),
             ((8 * g * x) / p + (9 * g**2 * x**2) - (8 * g * x**2)/p - (12 * g**2 * x**3)),
             (4/p + (6 * g * x) - (8 * x)/p - (18 * g * x**2) + (4 * x**2)/p + (12 * g * x**3)),
             -6 * x + 9 * x**2 - 4*x**3]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    return max(result)


def wlc(distances, length, p, method=marko_siggia, k=0, residues_distance=None):
    if not residues_distance:
        residues_distance = default_parameters['residues_distance']
    if method == 'stretch-adjusted':
        return np.array([stretch_adjusted_wlc(d, length, p, k, residues_distance) for d in distances])
    else:
        return np.array([marko_siggia(d, length, p, k) for d in distances])


# calculating d(F)
def inverse_marko_siggia(f, p, k=0):
    if k == 0:
        coefs = [1, -(2.25 + f / p), (1.5 + 2 * f / p), -f / p]
    else:
        coefs = [1,
                 -(2.25 + f * (3*k + 1 / p)),
                 (3 * (k ** 2) + 2 * (k / p)) * f ** 2 + ((4.5 * k) + (2 / p)) * f + 1.5,
                 -f * (((k**3) + ((k**2) / p)) * (f ** 2) + (2.25 * (k**2) + 2 * (k / p)) * f + ((1.5 * k) + (1 / p)))]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    if k == 0:
        result = result[result < 1]
    return min(result)


def inverse_stretch_adjusted_wlc(f, p, k, residues_distance=None):
    if not residues_distance:
        residues_distance = default_parameters['residues_distance']
    g = k * residues_distance
    coefs = [(-4 + 12 * f * g - 12 * f**2 * g**2 + 4 * f**3 * g**3),
             (9 - 18 * f * g + 9 * f**2 * g**2 + (4 * f)/p - (8 * f**2 * g)/p + (4 * f**3 * g**2)/p),
             (-6 + 6 * f * g - (8 * f)/p + (8 * f**2 * g)/p),
             4 * f/p]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    return min(result)


def inverse_wlc(forces, p, method='marko_siggia', k=0, residues_distance=None):
    if not residues_distance:
        residues_distance = default_parameters['residues_distance']
    if method == 'stretch-adjusted':
        return np.array([inverse_stretch_adjusted_wlc(f, p, k, residues_distance) for f in forces])
    else:
        return np.array([inverse_marko_siggia(f, p, k) for f in forces])


def get_d_dna(p_dna, l_dna, k_dna, f_space, method='marko-siggia'):
    if l_dna > 0:
        return l_dna * inverse_wlc(f_space, p_dna, k=k_dna, method=method)
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


def plot_decomposed_histogram(position, data, bound, residue_distance):
    k = 0
    l_space = np.linspace(0, bound, 1001)
    # TODO add Cauchy distribution
    for index, row in data[['means', 'widths', 'heights', 'cauchy_means', 'cauchy_gammas', 'factor']].iterrows():
        mean, width, height, cauchy_means, cauchy_gammas, factor = tuple(row.to_numpy())
        residues = 1 + int(mean / residue_distance)
        label = "L= " + str(round(mean, 3)) + ' (' + str(residues) + ' AA)'
        # normal distribution
        # y_gauss = factor * gauss(l_space, height, mean, width)
        # position.plot(l_space, y_gauss, linestyle='--', linewidth=0.5, label=label, color=get_color(k, len(data)))
        # cauchy distribution
        y_cauchy = factor * cauchy.pdf(l_space, cauchy_means, cauchy_gammas)
        position.plot(l_space, y_cauchy, linestyle='--', linewidth=0.5, label=label, color=get_color(k, len(data)))
        k += 1
    return


def plot_trace_fits(position, coefficients, max_f, residue_distance, method='marko-siggia'):
    f_space = np.linspace(0.1, max_f)
    d_dna = get_d_dna(coefficients.get('p_dna', 0), coefficients.get('l_dna', 0),
                      coefficients.get('k_dna', None), f_space, method=method)
    x_prot = inverse_wlc(f_space, coefficients.get('p_prot', 0), k=coefficients.get('k_prot', None), method=method)
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
def find_state_boundaries(smooth_data, parameters, p, k, max_distance=0.3, method='marko-siggia'):
    begs = []
    ends = []
    last_end = 0
    for index, row in parameters.iterrows():
        l_prot = row['means']
        smooth_data['state_' + str(index)] = wlc(list(smooth_data['d']), l_prot, p, k=k, method=method)
        data_close = smooth_data[abs(smooth_data['F'] - smooth_data['state_' + str(index)]) <= max_distance]
        data_close = data_close[data_close['d'] > last_end]['d']
        begs.append(data_close.min())
        ends.append(data_close.max())
        if ends[-1] != np.NaN:
            last_end = ends[-1]
    return begs, ends, smooth_data


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


def get_max_contour_length(data, significance):
    if isinstance(data, np.ndarray):
        series = pd.Series(data)
    elif isinstance(data, pd.DataFrame):
        series = data['hist_values']
    else:
        series = data
    hist_values = series.astype(int).value_counts(normalize=True)
    hist_values = hist_values[hist_values > significance].sort_index()
    max_contour_length = hist_values.index[-1] + 20
    return max_contour_length


def decompose_histogram(hist_values, significance=None, states=None):
    if not significance:
        significance = default_parameters['significance']

    if not states:
        maximas, support = guess_states_number(hist_values, significance=significance)
        states = max(len(maximas), 1)
    else:
        support = [0, get_max_contour_length(hist_values, significance)]

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
    skewness = []
    factors = []
    products = []
    means = np.array(parameters['means'].values)
    middles = means[:-1] + 0.5 * np.diff(means)
    boundaries = [0] + list(middles) + [max(trimmed_data)]
    for k in range(len(states_names)):
        state = states_names[k]
        bounds = [boundaries[k], boundaries[k+1]]
        matching = states[(states['state'] == state) & (states['d'].between(bounds[0], bounds[1]))]['d'].to_numpy()
        factors.append(float(len(matching))/len(states))
        mu, gamma = cauchy.fit(matching)
        ensamble = cauchy.rvs(size=len(matching), loc=mu, scale=gamma)
        cauchy_means.append(mu)
        cauchy_gammas.append(gamma)
        skewness.append(skew(matching))
        pvalues.append(ks_2samp(ensamble, matching).pvalue)
    parameters['cauchy_means'] = np.array(cauchy_means)
    parameters['cauchy_gammas'] = np.array(cauchy_gammas)
    parameters['skewness'] = np.array(skewness)
    parameters['pvalues'] = np.array(pvalues)
    parameters['factor'] = np.array(factors)
    return parameters, support


def valid_parameters(parameters):
    bounds_template = {'p': (0, 100), 'k': (0, 0.1), 'l': (0, np.inf)}
    for key in parameters.keys():
        if parameters[key] < bounds_template[key[0]][0] or parameters[key] > bounds_template[key[0]][1]:
            return False
    return True


# fitting
def fit_skew(x, data, states, known, unknown_keys, method='marko-siggia'):
    unknown = {unknown_keys[k]: x[k] for k in range(len(unknown_keys))}
    parameters = {**known, **unknown}

    # validation
    if not valid_parameters(parameters):
        return 999

    # data calculation
    d_dna = get_d_dna(parameters.get('p_dna'), parameters.get('l_dna'), parameters.get('k_dna'), (data['F']),
                      method=method)
    x_prot = inverse_wlc(data['F'], parameters.get('p_prot'), k=parameters.get('k_prot'), method=method)
    hist_values = (data['d'] - d_dna)/x_prot

    # parameter calculation
    parameters, support = decompose_histogram(hist_values, states=states)
    skewness = parameters['skewness'].abs().mean()
    return skewness


def fit_coefficients(data, parameters):
    # templates
    coefficients = ['p_prot', 'k_prot', 'p_dna', 'k_dna', 'l_dna']
    linker_known = {'p_dna': 0, 'k_dna': 0, 'l_dna': 0}

    # setting known values
    states = parameters['states']
    method = parameters['method']
    known = {c: parameters[c] for c in coefficients if parameters[c] >= 0}
    if not parameters['linker'] and not bool({'p_dna', 'l_dna', 'k_dna'} & set(known.keys())):
        known = {**known, **linker_known}

    # setting unknown values with initial guesses and bounds
    unknown = {c: parameters['initial_guess'][c] for c in coefficients if c not in known.keys()}
    if len(unknown.keys()) == 0:
        return known
    x0 = np.array(list(unknown.values()))

    x_opt = minimize(fit_skew, x0=x0, args=(data, states, known, list(unknown.keys()), method), method='Nelder-Mead')
    fitted = {list(unknown.keys())[k]: x_opt.x[k] for k in range(len(list(unknown.keys())))}
    result = {**known, **fitted}
    return result
