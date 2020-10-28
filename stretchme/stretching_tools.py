""" The miscellanea functions for stretching package"""

import pandas as pd
import numpy as np
import logging
import matplotlib.colors as mcolors
from os import path
from io import StringIO
from scipy.optimize import curve_fit, minimize, fsolve
from scipy.special import erf
from scipy.stats import skew
from scipy.signal import argrelextrema
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
    if columns and all([isinstance(c, int) for c in columns]):
        return input_data.iloc[:, columns]
    elif columns and not all([isinstance(c, int) for c in columns]):
        return input_data[columns]
    elif cases:
        allowed = [str(_) for _ in cases]
        colnames = [name for name in list(input_data) if name.strip('dF_') in allowed]
        return input_data[colnames]
    else:
        return input_data


def read_from_file(input_data, cases, columns, parameters, name, debug):
    print(input_data)
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
def single_gaussian(x, height, mean, width):
    if height < 0 or mean < 0:
        return np.zeros(len(x))
    return height * np.exp(-(x - mean)**2/(2*width**2))


def multiple_gaussian(x, *args):
    """ The function expects 3*states parameters (height1, center1, width1, height2, center2, width2, ..."""
    result = args[-1] * np.ones(len(x))
    for k in range(0, len(args)-1, 3):
        height, mean, width = args[k], args[k + 1], args[k + 2]
        result += single_gaussian(x, height, mean, width)
    return result


def integrate_gauss(force, mean, width, height):
    # return width * np.sqrt(2*np.pi) * (1 - 0.5 * erf((force - mean)/(width * np.sqrt(2))))
    return 0.5 * (1-erf((force - mean)/(np.sqrt(2) * width)))

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


def func(f, p, k, length, d):
    first_part = 0.25 / ((1 - d / length + f * k) ** 2) - 0.25 + d / length - f * k
    return p * (first_part - 0.8 * np.power(d / length - f * k, 2.15)) - f


def improved_marko_siggia(d, length, p, k=0):
    if k == 0 and d > 0.99 * length:
        return 999
    elif k == 0:
        return p * (0.25 / ((1 - d / length) ** 2) - 0.25 + d / length - 0.8 * np.power(d / length, 2.15))
    else:
        return fsolve(func, marko_siggia(d, length, p, k), args=(p, k, length, d))


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
    result = result[result >= 0]
    return max(result)


def wlc(distances, length, p, method='marko_siggia', k=0, residues_distance=None):
    if not residues_distance:
        residues_distance = default_parameters['residues_distance']
    if method == 'stretch-adjusted':
        return np.array([stretch_adjusted_wlc(d, length, p, k, residues_distance) for d in distances])
    if method == 'improved_marko_siggia':
        return np.array([improved_marko_siggia(d, length, p, k) for d in distances])
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


def inverse_improved_marko_siggia(f, p, k=0):
    exponent = np.power(900 * p / f, 0.25)
    if k == 0:
        return 4 / 3 - 4 / (3 * np.sqrt(f / p + 1)) - \
               10 * np.exp(exponent) / (np.sqrt(f / p) * (np.exp(exponent) - 1)**2) + \
               np.power(f / p, 1.62) / (3.55 + 3.8 * np.power(f / p, 2.2))
    else:
        return 4 / 3 - 4 / (3 * np.sqrt(f / p + 1)) - \
               10 * np.exp(exponent) / (np.sqrt(f / p) * (np.exp(exponent) - 1)**2) + \
               np.power(f / p, 1.62) / (3.55 + 3.8 * np.power(f / p, 2.2)) + f * k


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
    if method == 'improved_marko_siggia':
        return np.array([inverse_improved_marko_siggia(f, p, k) for f in forces])
    else:
        return np.array([inverse_marko_siggia(f, p, k) for f in forces])


def get_d_dna(p_dna, l_dna, k_dna, f_space, method='marko-siggia'):
    if l_dna > 0:
        return l_dna * inverse_wlc(f_space, p_dna, k=k_dna, method=method)
    else:
        return np.zeros(len(f_space))


def get_force_load(force, p_linker, l_linker, speed, k_spring):
    # speed = parameters['speed']
    # l_dna, p_dna, k_dna = parameters['l_dna'], parameters['p_dna'], parameters['k_dna']
    # k_spring = parameters['spring_constant']
    if float(k_spring) == 0.0:
        raise NotImplementedError("Cannot deal with no-spring case right now")
    dna_part = (2 * (l_linker/p_linker) * (1 + force/p_linker)) / (3 + 5 * force/p_linker + 8 * (force/p_linker)**(5/2))
    factor = dna_part + 1/k_spring
    return speed/factor


def dhs_feat_cusp(force, x, t0, g):
    if 1 - 0.5 * force.max() * x/g < 0 or t0 < 0:
        return np.array([999 for _ in range(len(force))])
    return np.log(t0) - x*force - np.log(1 - 0.5 * force * x/g) + ((0.5*force*x)**2) / g


def dhs_feat_linear_cubic(force, x, t0, g):
    return t0 / (1 - 2 * x/g * force/3)**(1/2) * np.exp(-g*(1-(1-2 * x/g * force/3)**(3/2)))


def dhs_feat_bell(force, x, t0):
    if t0 < 0:
        return np.array([999 for _ in range(len(force))])
    return np.log(t0) - x*force


def dhs_feat(data, init_x):
    coefficients = {}
    init_lifetime = data['lifetime'].head(1).values[0]

    # v = 1
    p0 = (init_x, init_lifetime)
    popt, pcov = curve_fit(dhs_feat_bell, data['forces'], np.log(data['lifetime']), p0=p0)
    coefficients['bell'] = {'x': popt[0], 't0': popt[1]}   #'covariance': pcov}

    # v = 1/2
    p0 = (coefficients['bell']['x'], coefficients['bell']['t0'], coefficients['bell']['x']*data['forces'].max())
    try:
        popt, pcov = curve_fit(dhs_feat_cusp, data['forces'], np.log(data['lifetime']), p0=p0)
        result = {'x': popt[0], 't0': popt[1], 'g': popt[2]}   #, 'covariance': pcov}
    except RuntimeError:
        result = None
    coefficients['cusp'] = result

    # v = 2/3
    p0 = (coefficients['bell']['x'], coefficients['bell']['t0'], coefficients['bell']['x'] * data['forces'].max())
    try:
        popt, pcov = curve_fit(dhs_feat_linear_cubic, data['forces'], data['lifetime'], p0=p0)
        result = {'x': popt[0], 't0': popt[1], 'g': popt[2]}  # , 'covariance': pcov}
    except RuntimeError:
        result = None
    coefficients['linear_cubic'] = result
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
    print('data for histogram decomposition:')
    print(data)
    data = data[data['widths'] < 10]
    for index, row in data[['means', 'widths', 'heights']].iterrows():
        mean, width, height = tuple(row.to_numpy())
        residues = 1 + int(mean / residue_distance)
        label = "L= " + str(round(mean, 3)) + ' (' + str(residues) + ' AA)'
        y_gauss = single_gaussian(l_space, height, mean, width)
        position.plot(l_space, y_gauss, linestyle='--', linewidth=0.5, label=label, color=get_color(k, len(data)))
        k += 1
    return


def plot_trace_fits(position, coefficients, max_f, residue_distance, method='marko-siggia'):
    f_space = np.linspace(0.1, max_f)
    d_dna = get_d_dna(coefficients.get('p_dna', 0), coefficients.get('l_dna', 0),
                      coefficients.get('k_dna', None), f_space, method=method)
    x_prot = inverse_wlc(f_space, coefficients.get('p_prot', 0), k=coefficients.get('k_prot', None), method=method)
    k = 0
    parameters = coefficients['l_prot'][coefficients['l_prot']['widths'] < 10]
    for index, row in parameters.iterrows():
        l_prot = row['means']
        residues = 1 + int(l_prot / residue_distance)
        d_prot = l_prot * x_prot
        d_plot = d_prot + d_dna
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
    return begs, ends


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


def transform_guesses(guesses, background_level, min_mean, max_mean):
    p0 = []
    bound_min = tuple([item for k in range(len(guesses)) for item in [0, min_mean, 0]] + [0])
    bound_max = tuple([item for k in range(len(guesses)) for item in [np.inf, max_mean, np.inf]] + [np.inf])
    for ind, row in guesses.iterrows():
        p0 += list(row.values)
    p0.append(background_level)
    return tuple(p0), (bound_min, bound_max)


def transform_fit_results(results):
    parameters = pd.DataFrame({'heights': np.array([results[k] for k in range(0, len(results)-1, 3)]),
                               'means': np.array([results[k+1] for k in range(0, len(results)-1, 3)]),
                               'widths': np.array([abs(results[k+2]) for k in range(0, len(results)-1, 3)])})
    parameters = parameters.sort_values(by=['means'])
    return parameters


def decompose_histogram(hist_values, significance=None, states=None, bandwidth=None, max_value=None,
                        background_level=0, init_means=None):
    if not significance:
        significance = default_parameters['significance']

    if not bandwidth:
        bandwidth = default_parameters['bandwidth']

    # finding the density from histogram
    hist_values = hist_values[hist_values > 0]
    if max_value:
        hist_values = hist_values[hist_values < max_value]
    x = np.expand_dims(hist_values, 1)
    kde = KernelDensity(bandwidth=bandwidth).fit(x)
    estimator = np.linspace(min(hist_values), max(hist_values), 1001)
    kde_est = np.exp(kde.score_samples(estimator.reshape(-1, 1)))
    kde_est[kde_est < background_level] = background_level
    significant = [estimator[_] for _ in range(len(estimator)) if kde_est[_] > significance]
    support = [min(significant), max(significant)]
    trimmed_data = hist_values[(hist_values >= support[0]) & (hist_values <= support[1])]

    # finding the local maxima of the histogram = number of states and their heights and widths
    if init_means:
        means = np.array(init_means)
    else:
        means = np.array([estimator[_] for _ in argrelextrema(kde_est, np.greater)[0] if kde_est[_] > significance])
    if states:
        missing = max(states - len(means), 0)
        if missing > 0:
            intervals = missing + 1
            beg, end = min(means), max(means)
            additional = np.array([beg * (intervals - i)/intervals + end * i/intervals for i in range(1, intervals)])
            means = np.append(means, additional)
    heights = np.exp(kde.score_samples(means.reshape(-1, 1)))
    guesses = pd.DataFrame({'heights': heights, 'means': means, 'widths': np.ones(len(means))})
    guesses = guesses.sort_values(by=['heights'], ascending=False)
    if guesses['heights'].max() > 0.01:
        guesses = guesses[guesses['heights'] > 0.0075]

    # for small number of points:
    if len(hist_values) <= 5:
        states = 1

    guesses = guesses.head(states)
    # print(guesses)

    # fitting
    p0, bounds = transform_guesses(guesses, background_level, min(hist_values), max(hist_values))
    # print('bounds', bounds)
    popt = tuple()
    while len(popt) < len(p0):
        try:
            popt, pcov = curve_fit(multiple_gaussian, estimator, kde_est, p0=p0, bounds=bounds)
        except RuntimeError:
            p0 = tuple(list(p0[:-4]) + [p0[-1]])
            bounds = (tuple(list(bounds[0])[:-4] + [list(bounds[0])[-1]]),
                      tuple(list(bounds[1])[:-4] + [list(bounds[1])[-1]]))
            print("I reduced the number of states from expected (" + str(int((len(p0) + 3)/3)) + ") to " +
                  str(int(len(p0)/3)))
    if len(popt) <= 1:
        mean = sum(hist_values)/len(hist_values)
        width = 0.5 * (max(hist_values) - min(hist_values))
        popt = (1, mean, width)
        print("Putting the parameters " + str(popt))
    parameters = transform_fit_results(popt)

    #calculating the probability the value corresponds to a given state

    states = pd.DataFrame({'d': trimmed_data})
    states_names = []

    for ind, row in parameters.iterrows():
      states_names.append('state_' + str(ind))
      mean, width, height = row[['means', 'widths', 'heights']].values
      states[states_names[-1]] = single_gaussian(trimmed_data, height, mean, width)


    #assigning the best matching state for each value
    states['state'] = states[states_names].idxmax(axis=1)

    # searching for the begin and end of each state
    begs, ends, skewness = [], [], []
    for k in range(len(states_names)):
        state = states_names[k]
        matching = states[(states['state'] == state)]['d']
        beg = matching.min()
        if k > 0:
            beg = max(beg, ends[-1])
        begs.append(beg)
        ends.append(matching.max())
        skewness.append(skew(matching))
    print(begs, ends)
    parameters['begs'], parameters['ends'], parameters['skewness'] = begs, ends, skewness
    parameters = parameters[parameters['heights'] > 0.001]
    return parameters, support


def valid_parameters(parameters, max_prot, max_dna=np.inf, inverse=False):
    # bounds_template = {'p': (0.1, 100), 'k': (0, 0.1), 'l': (0, max_dna), 'L': (0, max_prot)}
    bounds_template = {'p_dna': (0.13, 0.21), 'p_prot': (0.01, 1.5), #'p_prot': (5.66, 5.98),     # in T=298 p_dna should be > 0.13
                       'k_dna': (0, 0.1), 'k_prot': (0, 0.1),
                       'l_dna': (330, max_dna), 'l_prot_first': (370, 500), 'l_prot_last': (370, 500)}    # 'l_dna': (330, max_dna), , 'l_prot': (0, max_prot)
    for key in parameters.keys():
        # if parameters['l_dna'] == 0 and 'dna' in key:
        #     continue
        # if parameters[key] < bounds_template[key[0]][0] or parameters[key] > bounds_template[key[0]][1]:
        if parameters[key] < bounds_template[key][0] or parameters[key] > bounds_template[key][1]:
            return False
    # P_L prot ~ 0.7, P_L DNA < 30 to (kbT/P_Lprot = p_prot) / (p_dna = kbT/P_L DNA) < 30/0.7 < 50
    if parameters['p_dna'] > 0 and parameters['p_prot'] / parameters['p_dna'] > 50:
            return False
    return True


# fitting
def find_range(data, data_smooth, first_range=None, last_range=None):
    minimas_indices = argrelextrema(data_smooth['F'].to_numpy(), np.less)[0]
    maximas_indices = argrelextrema(data_smooth['F'].to_numpy(), np.greater)[0]
    extremas_indices = [item for sublist in zip(minimas_indices, maximas_indices) for item in sublist] + \
                       [len(data_smooth) - 1]
    extremas = data_smooth.loc[extremas_indices, ['d', 'F']].values
    for k in range(1, len(extremas)):
        if extremas[k][0] - extremas[k-1][0] > 3 and extremas[k-1][1] - extremas[k][1] > 0.5:
            first_range = (data['d'].min(), extremas[k-1][0])
            break
    for l in range(len(extremas)-1, 0, -1):
        if extremas[l][0] - extremas[l-1][0] > 3 and extremas[l][1] - extremas[l-1][1] > 5:
            last_range = (extremas[l-1][0], extremas[l][0])
            break
    return (first_range, last_range)


def fit_last(x, data, last_range, max_length, known, unknown_keys, method='marko-siggia', residues_distance=None):
    print(x)
    unknown = {unknown_keys[k]: x[k] for k in range(len(unknown_keys))}
    parameters = {**known, **unknown}

    # validation
    if not valid_parameters(parameters, max_length, 370):
        # print("invalid")
        return np.inf

    fit_data = data[data['d'].between(last_range[0], last_range[1])]

    if parameters['l_dna'] > 0:
        d_dna = get_d_dna(parameters['p_dna'], parameters['l_dna'], parameters['k_dna'], fit_data['F'].to_numpy(),
                          method=method)
        fit_data['d'] = fit_data['d'] - d_dna
        fit_data = fit_data[fit_data['d'] > 0]
        fit_data = fit_data[fit_data['d'] < parameters['l_prot']]


    fit_f = wlc(fit_data['d'], parameters['l_prot'], parameters['p_prot'], method=method, k=parameters['k_prot'],
                residues_distance=residues_distance)
    return np.linalg.norm(fit_f - fit_data['F'].to_numpy())


def fit_parts(x, data, first_range, last_range, max_length, known, unknown_keys, method='marko-siggia',
              residues_distance=None):
    print(x)

    unknown = {unknown_keys[k]: x[k] for k in range(len(unknown_keys))}
    parameters = {**known, **unknown}

    # validation
    if not valid_parameters(parameters, max_length, 370):
        # print("invalid")
        return np.inf

    fit_first = data[data['d'].between(first_range[0], first_range[1])]
    fit_last = data[data['d'].between(last_range[0], last_range[1])]

    fit_value_first = wlc(fit_first['d'], parameters['l_first'], parameters['p_tot'], method=method,
                          k=parameters['k_tot'], residues_distance=residues_distance)
    fit_value_last = wlc(fit_last['d'], parameters['l_tot'], parameters['p_tot'], method=method,
                         k=parameters['k_tot'], residues_distance=residues_distance)

    error = np.linalg.norm(fit_value_first - fit_first['F'].to_numpy()) + \
            np.linalg.norm(fit_value_last - fit_last['F'].to_numpy())

    return error


def fit_coefficients(data, data_smooth, parameters, experiment, debug=False):
    # templates
    coefficients = ['p_tot', 'l_tot', 'k_tot']  #['p_prot', 'k_prot', 'p_dna', 'k_dna', 'l_dna']
    # linker_known = {'p_dna': 0, 'k_dna': 0, 'l_dna': 0}

    # setting known values
    method = parameters['method']
    residues = parameters['residues']
    residues_distance = parameters['residues_distance']
    first_range = parameters['first_range']
    last_range = parameters['last_range']
    if not first_range or not last_range:
        first_range, last_range = find_range(data, data_smooth, first_range, last_range)
    print(first_range, last_range)
    if debug:
        logger = set_logger(experiment)
        logger.info("Ranges:\n" + str(first_range) + '\t' + str(last_range))
        close_logs(logger)
    known = {c: parameters.get(c, 0) for c in coefficients if parameters.get(c, 0) > 0}
    max_length = residues_distance * (residues - 1) + 370
    print(residues, max_length)

    # if not parameters['linker'] and not bool({'p_dna', 'l_dna', 'k_dna'} & set(known.keys())):
    #     known = {**known, **linker_known}

    # setting unknown values with initial guesses and bounds
    unknown = {c: parameters['initial_guess'][c] for c in coefficients if c not in known.keys()}
    if len(unknown.keys()) == 0:
        return known
    guess_prot_length = residues_distance * (residues - 1 - parameters['missing']) + 340    #parameters['l_dna']
    unknown['l_tot'] = guess_prot_length
    unknown['l_first'] = guess_prot_length - 50
    print(unknown)

    x0 = np.array(list(unknown.values()))
    x_opt = minimize(fit_parts, x0=x0, args=(data, first_range, last_range, max_length, known, list(unknown.keys()),
                                             method, None), method='Nelder-Mead')

    fitted = {list(unknown.keys())[k]: x_opt.x[k] for k in range(len(list(unknown.keys())))}
    result = {**known, **fitted}
    return result


def extract_dna(data, parameters):
    return


def fit_coefficients2(data, l_tot, method='marko-siggia'):
    return