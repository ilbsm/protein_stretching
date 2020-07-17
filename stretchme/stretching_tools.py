""" The miscellanea functions for stretching package"""

import pandas as pd
import numpy as np
import logging
import matplotlib.colors as mcolors
from os import path
from io import StringIO
from scipy.optimize import curve_fit, minimize
from scipy.special import erf
from scipy.stats import cauchy, ks_2samp
from scipy.signal import argrelextrema
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from colorsys import hsv_to_rgb
from .default_parameters import default_parameters


# preprocessing
def pack_parameters(filename, sheet_name=0, linker=None, unit='nm', speed=1, residues_distance=0.365, source='theory',
                    low_force_cutoff=0.1, plot_columns=4, initial_guess=None, separator=None):

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


def ewlc_np(d, p, k, length):
    return np.array([reduced_ewlc(point/length, p, k) for point in d])


def wlc_np(d, p, length):
    return np.array([reduced_wlc(point/length, p) for point in d])


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


def plot_decomposed_histogram(position, data, bound, residue_distance):
    k = 0
    l_space = np.linspace(0, bound, 1001)
    # TODO add Cauchy distribution
    for index, row in data[['means', 'widths', 'heights']].iterrows():
        mean, width, height = tuple(row.to_numpy())
        y_plot = gauss(l_space, height, mean, width)
        residues = 1 + int(mean / residue_distance)
        label = "L= " + str(round(mean, 3)) + ' (' + str(residues) + ' AA)'
        position.plot(l_space, y_plot, linestyle='--', label=label, color=get_color(k, len(data)))
        k += 1
    return


def plot_trace_fits(position, coefficients, max_f, residue_distance):
    f_space = np.linspace(0.1, max_f)
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


def decompose_histogram(hist_values, significance=0.01, states=None,):
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
    return parameters, support


# fitting
def fit_error(x, data, states, known, unknown_keys):
    unknown = {unknown_keys[k]: x[k] for k in range(len(unknown_keys))}
    parameters = {**known, **unknown}
    d_dna = parameters.get('l_dna') * invert_wlc_np(data['F'], parameters.get('p_dna'), parameters.get('k_dna'))
    x_prot = invert_wlc_np(data['F'], parameters.get('p_prot'), parameters.get('k_prot'))
    hist_values = (data['d'] - d_dna)/x_prot
    try:
        parameters, support = decompose_histogram(hist_values, states=states)
        return parameters['cauchy_gammas'].min() - parameters['pvalues'].max()
    except:
        return 999


def fit_coefficients(data, parameters):
    # templates
    coefficients = ['p_prot', 'k_prot', 'p_dna', 'k_dna', 'l_dna']
    linker_known = {'p_dna': 0, 'k_dna': 0, 'l_dna': 0}
    bounds_template = {'p': (0.1, 100), 'k': (0, 1), 'l': (0, np.inf)}

    # setting known values
    known = {c: parameters[c] for c in coefficients if parameters[c] >= 0}
    if not parameters['linker'] and not bool({'p_dna', 'l_dna', 'k_dna'} & set(known.keys())):
        known = {**known, **linker_known}

    # setting unknown values with initial guesses and bounds
    unknown = {c: parameters['initial_guess'][c] for c in coefficients if c not in known.keys()}
    if len(unknown.keys()) == 0:
        return known
    x0 = np.array(list(unknown.values()))
    bounds = tuple([bounds_template[key[0]] for key in unknown.keys()])

    x_opt = minimize(fit_error, x0=x0, args=(data, parameters['states'], known, list(unknown.keys())), bounds=bounds)
    fitted = {list(unknown.keys())[k]: x_opt.x[k] for k in range(len(list(unknown.keys())))}
    result = {**known, **fitted}
    return result
