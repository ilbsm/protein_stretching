""" The miscelenea functions for stretching package"""

import pandas as pd
import numpy as np
from .default_parameters import default_parameters
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
from colorsys import hsv_to_rgb

colors = [mcolors.CSS4_COLORS['red'],
          mcolors.CSS4_COLORS['green'],
          mcolors.CSS4_COLORS['blue'],
          mcolors.CSS4_COLORS['yellow'],
          mcolors.CSS4_COLORS['cyan'],
          mcolors.CSS4_COLORS['orange'],
          mcolors.CSS4_COLORS['purple'],
          mcolors.CSS4_COLORS['lime'],
          mcolors.CSS4_COLORS['magenta']]


def find_derivative(x, y):
    y_diff = np.diff(np.array(y)) / np.diff(np.array(x))
    x_diff = np.array(x[:-1]) + np.diff(np.array(x)) / 2
    return running_average(x_diff, y_diff)


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


def cluster_coefficients(coefficients, maxgap=15, minnumber=4, minspan=-1):
    coefficients.sort()
    clusters = [[coefficients[0]]]
    for x in coefficients[1:]:
        if abs(x - clusters[-1][-1]) <= maxgap:
            clusters[-1].append(x)
        else:
            clusters.append([x])
    singles = [k for k in range(len(clusters)) if len(clusters[k]) <= minnumber or
               (max(clusters[k]) - min(clusters[k])) < minspan]
    if singles:
        clusters.append([loc for k in singles for loc in clusters[k]])
        for k in list(reversed(singles)):
            clusters.pop(k)
    return clusters


def running_average(x, y, window=None):
    if not window:
        window = max(int(len(x) / 100), 8)
    x_smooth = np.convolve(x, np.ones((window,)) / window, mode='valid')
    y_smooth = np.convolve(y, np.ones((window,)) / window, mode='valid')
    return x_smooth, y_smooth


def invert_wlc(f, p, k=None):
    if not k:
        coefs = [1, -(2.25 + f / p), (1.5 + 2 * f / p), -f / p]
    else:
        coefs = [1,
                 -(2.25 + f * (3 / k + 1 / p)),
                 (3 / k ** 2 + 2 / (k * p)) * f ** 2 + (4.5 / k + 2 / p) * f + 1.5,
                 -f * ((1 / k ** 3 + 1 / (p * k ** 2)) * f ** 2 + (2.25 / k ** 2 + 2 / (k * p)) * f + (
                             1.5 / k + 1 / p))]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    if not k:
        result = result[result < 1]
    return min(result)


def wlc(d, l, p, k):
    if not k:
        if d < 0.99 * l:
            return p * (0.25 / ((1 - d / l) ** 2) - 0.25 + d / l)
        else:
            return 999
    else:
        x = d/l
        coefs = [-(1/k**3) - 1/(p*k**2),
                 -(2.25/k**2) - 2/(k*p) + x * (3/k**2 + 2/(k*p)),
                 -(1.5/k) - 1/p + x * (4.5/k + 2/p) + x**2 * (-(3/k) - 1/p),
                 1.5 * x - 2.25 * x**2 + x**3]
        result = np.roots(coefs)
        result = np.real(result[np.isreal(result)])
        result = result[result > 0]
        return max(result)


def find_area(n, bins):
    width = bins[1] - bins[0]
    return sum(n) * width


def implicite_elastic_wlc(data, l, k, p):
    d = np.array(data['d'])
    f = np.array(data['F'])
    x = d / l
    coefficients = [1,
                    -f * (3 / k + 1 / p) - 9 / 4,
                    f ** 2 * (3 / (k ** 2) + 2 / (k * p)) + f * (9 / (2 * k) + 2 / p) + 3 / 2,
                    -f ** 3 * (1 / (k ** 3) + 1 / (p * k ** 2)) - f ** 2 * (9 / (4 * k ** 2) + 2 / (k * p)) - f * (
                                3 / (2 * k) + 1 / p)]
    return x ** 3 * coefficients[0] + x ** 2 * coefficients[1] + x * coefficients[2] + coefficients[3]


def implicite_elastic_wlc_amplitude(data, l, k, p):
    return np.abs(implicite_elastic_wlc(data, l, k, p))


def minimize_kp(df, length, init_k, init_p):
    ydata = np.zeros(len(df))
    popt, pcov = curve_fit(implicite_elastic_wlc_amplitude, df, ydata, bounds=(0, np.inf), p0=(length, init_k, init_p))
    return popt, pcov


def decompose_histogram(data, boundary_value):
    # finding the number of components
    x = np.expand_dims(data, 1)
    kde = KernelDensity().fit(x)
    estimator = np.linspace(min(data), max(data), 1001)
    kde_est = np.exp(kde.score_samples(estimator.reshape(-1, 1)))
    maximas = [estimator[_] for _ in argrelextrema(kde_est, np.greater)[0] if kde_est[_] > boundary_value]

    # finding the range
    # TODO clean it up
    min_boundary = max(
        [estimator[_] for _ in range(len(estimator)) if kde_est[_] < boundary_value and estimator[_] < maximas[0]])
    max_boundary = min(
        [estimator[_] for _ in range(len(estimator)) if kde_est[_] < boundary_value and estimator[_] > maximas[-1]])

    # fitting Gaussians
    x = np.expand_dims(data[(data > min_boundary) & (data < max_boundary)], 1)
    gmm = GaussianMixture(n_components=len(maximas))  # gmm for two components
    gmm.fit(x)
    # finding parameters
    parameters = pd.DataFrame({'means': np.array([x[0] for x in gmm.means_]),
                               'widths': np.array([x[0][0] for x in gmm.covariances_]),
                               'heights': np.array([np.exp(gmm.score_samples(np.array(u).reshape(-1, 1)))[0]
                                                    for u in [x[0] for x in gmm.means_]])})
    parameters = parameters.sort_values(by=['means'])
    boundaries = [min_boundary, max_boundary]
    return parameters, boundaries


def gauss(x, height, mean, width):
    return height * np.exp(-((x - mean) ** 2) / (2 * width))


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


def get_d_dna(p_dna, l_dna, k_dna, f_space):
    if p_dna > 0:
        return l_dna * np.array([invert_wlc(f, p_dna, k_dna) for f in f_space])
    else:
        return np.zeros(len(f_space))


def get_gray_shade(k, max_value):
    values = np.linspace(0, 0.75, max_value)
    return hsv_to_rgb(0, 0, values[k])


def get_color(k, max_value):
    values = np.linspace(0, 0.66, max_value)
    return hsv_to_rgb(values[k], 1, 1)


