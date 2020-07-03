""" The miscelenea functions for stretching package"""

import pandas as pd
import numpy as np
from .default_parameters import default_parameters
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from scipy.special import erf
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
from colorsys import hsv_to_rgb
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d


# preprocessing
def pack_parameters(filename, sheet_name=0, residues=None, distance=None, linker=None, source=None, unit=None,
                    speed=None, residues_distance=0.365, minimal_stretch_distance=10, high_force_cutoff=None,
                    low_force_cutoff=None, max_rupture_force=None, max_cluster_gap=15, plot_columns=4,
                    initial_guess=None,
                    separator=None):
    """ Filtering and packing the parameters for a clearer view."""
    # TODO add filters of the input parameters
    parameters = {
        'filename': filename,
        'sheet_name': sheet_name,
        'residues': residues,
        'distance': distance,
        'linker': linker,
        'source': source,
        'unit': unit,
        'speed': speed,
        'residues_distance': residues_distance,
        'minimal_stretch_distance': minimal_stretch_distance,
        'max_cluster_gap': max_cluster_gap,
        'plot_columns': plot_columns,
        'separator': separator
    }
    if high_force_cutoff:
        parameters['high_force_cutoff'] = high_force_cutoff
    else:
        parameters['high_force_cutoff'] = default_parameters['high_force_cutoff'][parameters['source']]
    if low_force_cutoff:
        parameters['low_force_cutoff'] = low_force_cutoff
    else:
        parameters['low_force_cutoff'] = default_parameters['low_force_cutoff'][parameters['source']]
    if max_rupture_force:
        parameters['max_rupture_force'] = max_rupture_force
    else:
        parameters['max_rupture_force'] = default_parameters['max_rupture_force'][parameters['source']]
    if high_force_cutoff:
        parameters['initial_guess'] = initial_guess
    else:
        parameters['initial_guess'] = default_parameters['initial_guess']
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


def invert_wlc_np(forces, p, k=None):
    return np.array([invert_wlc(f, p, k) for f in forces])


def invert_wlc(force, p, k=None):
    if not k:
        coefs = [1, -(2.25 + force / p), (1.5 + 2 * force / p), -force / p]
    else:
        coefs = [1,
                 -(2.25 + force * (3 / k + 1 / p)),
                 (3 / k ** 2 + 2 / (k * p)) * force ** 2 + (4.5 / k + 2 / p) * force + 1.5,
                 -force * ((1 / k ** 3 + 1 / (p * k ** 2)) * force ** 2 + (2.25 / k ** 2 + 2 / (k * p)) * force + (
                             1.5 / k + 1 / p))]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    if not k:
        result = result[result < 1]
    return min(result)


def wlc(d, length, p, k):
    if not k:
        if d < 0.99 * length:
            return p * (0.25 / ((1 - d / length) ** 2) - 0.25 + d / length)
        else:
            return 999
    else:
        x = d/length
        return reduced_wlc(x, p, k)


def reduced_wlc(x, p, k):
    coefs = [-(1/k**3) - 1/(p * k**2),
             -(2.25/k**2) - 2/(k * p) + x * (3/k**2 + 2/(k * p)),
             -(1.5/k) - 1/p + x * (4.5/k + 2/p) + x**2 * (-(3/k) - 1/p),
             1.5 * x - 2.25 * x**2 + x**3]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    return max(result)


def reduced_wlc_np(x, p, k):
    return np.array([reduced_wlc(point, p, k) for point in x])


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
    for index, row in data[['means', 'widths', 'heights']].iterrows():
        mean, width, height = tuple(row.to_numpy())
        y_plot = gauss(l_space, height, mean, width)
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
    parameters = pd.DataFrame({'means': np.array([x[0] for x in gmm.means_]),
                               'widths': np.array([x[0][0] for x in gmm.covariances_]),
                               'heights': np.array([np.exp(gmm.score_samples(np.array(u).reshape(-1, 1)))[0]
                                                    for u in [x[0] for x in gmm.means_]])})
    parameters = parameters.sort_values(by=['means'])
    boundaries = [min_boundary, max_boundary]
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


def minimize_pk(data, data_smoothed, p, k, significance=0.01, max_distance=0.01):
    hist_values = data['d']/np.array([invert_wlc(f, p, k) for f in data['F']])
    parameters, boundaries = decompose_histogram(np.array(hist_values), significance)
    smoothed = pd.DataFrame({'d': data_smoothed['d'], 'F': data_smoothed['F']})
    for index, row in parameters.iterrows():
        l_prot = row['means']
        smoothed['state_' + str(index)] = np.array([wlc(d, l_prot, p, k) for d in list(data_smoothed['d'])])
    last_end = 0
    for ind, row in parameters.iterrows():
        l_prot = row['means']
        data_close = smoothed[abs(smoothed['F'] - smoothed['state_' + str(ind)]) <= max_distance]
        data_close = data_close[data_close['d'] > last_end]['d']
        beg = data_close.min()
        end = data_close.max()
        last_end = end
        if l_prot == parameters['means'].max():
            fit_data = pd.DataFrame({'d': data[data['d'].between(beg, end)]['d'],
                                     'F': data[data['d'].between(beg, end)]['F'],
                                     'x': data[data['d'].between(beg, end)]['d']/l_prot})
            print(beg, end, l_prot)
    popt, pcov = curve_fit(reduced_wlc_np, fit_data['x'], fit_data['F'], p0=(p, k))
    print(popt)
    print(pcov)

    return 0.7735670704545268, 200


def implicit_elastic_wlc(data, p, k):
    result = pd.DataFrame(columns=['0', '1', '2', '3', 'sum'])
    result['3'] = data['x'] ** 3
    result['2'] = (-1) * data['x'] ** 2 * (2.25 + data['F'] * (3 / k + 1 / p))
    result['1'] = data['x'] * (data['F'] ** 2 * (3 / (k**2) + 2 / (k * p)) + data['F'] * (4.5 / k + 2 / p) + 1.5)
    result['0'] = data['F'] * ((1 / (k**3) + 1 / (p * k**2)) * data['F'] ** 2 +
                               (2.25 / (k**2) + 2 / (k * p)) * data['F'] + (1.5 / k + 1 / p))
    result['sum'] = np.abs(result['3'] + result['2'] + result['1'] + result['0'])
    return result['sum']



