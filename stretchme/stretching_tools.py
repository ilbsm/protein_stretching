""" The miscelenea functions for stretching package"""

import pandas as pd
import numpy as np
from .default_parameters import default_parameters
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors


colors = [mcolors.CSS4_COLORS['red'],
          mcolors.CSS4_COLORS['green'],
          mcolors.CSS4_COLORS['blue'],
          mcolors.CSS4_COLORS['yellow'],
          mcolors.CSS4_COLORS['cyan'],
          mcolors.CSS4_COLORS['magenta'],
          mcolors.CSS4_COLORS['orange'],
          mcolors.CSS4_COLORS['purple'],
          mcolors.CSS4_COLORS['lime']]


def find_derivative(x, y):
    y_diff = np.diff(np.array(y))/np.diff(np.array(x))
    x_diff = np.array(x[:-1]) + np.diff(np.array(x))/2
    return running_average(x_diff, y_diff)


def prepare_headers(data):
    """ Preparing the headers for dataframes"""
    headers = list(data)
    result = {}
    number = -1
    for k in range(len(headers)):
        if isinstance(headers[k], int) and k % 2 == 0:
            number = headers[k]
            result[headers[k]] = 'd_' + str(number)
        elif 'd' in headers[k] and k % 2 == 0:
            number = int(headers[k].strip('d_'))
            result[headers[k]] = 'd_' + str(number)
        elif k % 2 == 0:
            number = int(k/2)
            result[headers[k]] = 'd_' + str(number)
        else:
            result[headers[k]] = 'F_' + str(number)
    return result


def pack_parameters(filename, sheet_name=0, residues=None, distance=None, linker=None, source=None, unit=None,
                    speed=None, residues_distance=0.365, minimal_stretch_distance=10, high_force_cutoff=None,
                    low_force_cutoff=None, max_rupture_force=None, max_cluster_gap=15, columns=4, initial_guess=None,
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
        'columns': columns,
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


def read_data(filename, cases, sheet_name, separator):
    if '.xls' in filename:
        data = pd.read_excel(filename, sheet_name=sheet_name)
    elif '.csv' in filename:
        data = pd.read_csv(filename, separator=separator)
    else:
        raise NotImplemented("Not implemented format of the input data. Reading xls and csv files now.")
    data.dropna(axis='columns', how='all', inplace=True)
    headers = prepare_headers(data)
    data.rename(columns=headers, inplace=True)
    if cases:
        to_remove = list(data)
        for case in cases:
            to_remove.pop(to_remove.index('d_' + str(case)))
            to_remove.pop(to_remove.index('F_' + str(case)))
        data.drop(to_remove, axis=1, inplace=True)
    return data


def cluster_coefficients(coefficients, maxgap=15, minnumber=4, minspan=-1):
    coefficients.sort()
    clusters = [[coefficients[0]]]
    for x in coefficients[1:]:
        if abs(x - clusters[-1][-1]) <= maxgap:
            clusters[-1].append(x)
        else:
            clusters.append([x])
    singles = [k for k in range(len(clusters)) if len(clusters[k]) <= minnumber or
               (max(clusters[k])-min(clusters[k])) < minspan]
    if singles:
        clusters.append([loc for k in singles for loc in clusters[k]])
        for k in list(reversed(singles)):
            clusters.pop(k)
    return clusters


def running_average(x, y, window=None):
    if not window:
        window = int((max(x)-min(x))/100)
    x_smooth = np.convolve(x, np.ones((window,))/window, mode='valid')
    y_smooth = np.convolve(y, np.ones((window,))/window, mode='valid')
    return x_smooth, y_smooth


def find_hist_ranges(data):
    # TODO sprawdzic, czy wszystko jest dobrze z indeksami i dlugosciami list
    n, bins, patches = plt.hist(data, bins=200, density=False)
    plt.close()

    xs = []
    for k in range(len(bins) - 1):
        xs.append((bins[k] + bins[k + 1]) / 2)

    x_smooth, n_smooth = running_average(xs, n, window=5)
    # plt.plot(x_smooth, n_smooth)

    diff_x = np.diff(x_smooth)
    diff_n = np.diff(n_smooth)
    n_deriv = diff_n / diff_x
    n_bis = np.diff(diff_n) / diff_x[1:]

    to_cluster = []
    for k in range(len(n_bis)):
        if abs(n_deriv[k]) < 30 and n_bis[k] > 0:                        # not dna n_deriv < 3
            to_cluster.append(x_smooth[k + 1])

    clusters = cluster_coefficients(to_cluster, minnumber=0, maxgap=5)  # not dna maxgap=20
    fit_ranges = [[min(bins)]]
    for k in range(1, len(clusters) - 1):
        cluster = clusters[k]
        pos = (max(cluster) + min(cluster)) / 2
        fit_ranges[-1].append(pos)
        fit_ranges.append([pos])
    fit_ranges[-1].append(max(bins))
    return fit_ranges


def invert_wlc(f, p, k=None):
    if not k:
        coefs = [1, -(2.25 + f / p), (1.5 + 2 * f / p), -f / p]
    else:
        coefs = [1,
                 -(2.25 + f * (3/k + 1/p)),
                 (3/k**2 + 2/(k*p)) * f**2 + (4.5/k + 2/p) * f + 1.5,
                 -f * ((1/k**3 + 1/(p*k**2)) * f**2 + (2.25/k**2 + 2/(k*p)) * f + (1.5/k + 1/p))]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result > 0]
    if not k:
        result = result[result < 1]
    return min(result)


def find_area(n, bins):
    width = bins[1]-bins[0]
    return sum(n)*width



