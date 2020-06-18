import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from symfit import Variable, Parameter, Fit, Model, exp
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


def merge_dicts(dict1, dict2):
    return {**dict1, **dict2}


def invert_wlc_np(f, p, k=None):
    if f == 0:
        return 0
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


def calc_forward_wlc(d_array, length, pprot):
    return pprot * (0.25/((1-d_array/length)**2) - 0.25 + d_array/length)


def running_average(x, y, window=None):
    if not window:
        window = int((max(x)-min(x))/100)
    x_smooth = np.convolve(x, np.ones((window,))/window, mode='valid')
    y_smooth = np.convolve(y, np.ones((window,))/window, mode='valid')
    return x_smooth, y_smooth


def find_derivative(x, y):
    y_diff = np.diff(np.array(y))/np.diff(np.array(x))
    x_diff = np.array(x[:-1]) + np.diff(np.array(x))/2
    return running_average(x_diff, y_diff)


def calc_inverse_wlc_dna(data, ldna=0, pdna=0, pprot=0, *args):
    if any([len(r) != len(args) for r in data]):
        raise SyntaxError("The size of the data does not match the number of protein lengths given"
                          " (" + str(len(args)) + ").")

    d_dna = np.zeros(len(data))
    d_prot = np.zeros(len(data))

    if ldna and pdna:
        d_dna = np.array([2 * ldna * invert_wlc_np(sum(f), pdna) for f in data])
    for k in range(len(args)):
        lprot = args[k]
        d_prot += np.array([lprot * invert_wlc_np(f[k], pprot) for f in data])

    return d_prot + d_dna


def calc_inverse_wlc_prot(data, pprot=0, *args):
    if any([len(r) != len(args) for r in data]):
        raise SyntaxError("The size of the data does not match the number of protein lengths given"
                          " (" + str(len(args)) + ").")

    d_prot = np.zeros(len(data))
    for k in range(len(args)):
        lprot = args[k]
        d_prot += np.array([lprot * invert_wlc_np(f[k], pprot) for f in data])

    return d_prot


def generate_vectors(dist, forces, ranges, low_force_cutoff):
    dist_vector = np.array([dist[_] for _ in range(len(dist)) if forces[_] > low_force_cutoff])
    forces_vector = np.zeros((len(forces), len(ranges)))
    for r in range(len(ranges)):
        current_range = ranges[r]
        for k in range(len(forces)):
            f = forces[k]
            d = dist[k]
            if f > low_force_cutoff and current_range[0] <= d <= current_range[1]:
                forces_vector[k][r] = f
    return dist_vector, forces_vector


def wlc_fit(ranges, dist, forces, linker, low_force_cutoff):
    dist_vector, forces_vector = generate_vectors(dist, forces, ranges, low_force_cutoff)
    if linker == 'dna':
        bounds_min = [0 for k in range(len(ranges) + 3)]
        bounds_max = [np.inf, 20, 50]
        p0 = [200, 0.3, 30]
        for k in range(len(ranges)):
            p0.append(100*(k+1))
            bounds_max.append(np.inf)
        popt, pcov = curve_fit(calc_inverse_wlc_dna, forces_vector, dist_vector, p0=tuple(p0),
                               bounds=(bounds_min, bounds_max))
    else:
        bounds_min = [0 for k in range(len(ranges) + 1)]
        bounds_max = [50]
        p0 = [0]
        for k in range(len(ranges)):
            p0.append(100*(k+1))
            bounds_max.append(np.inf)
        popt, pcov = curve_fit(calc_inverse_wlc_prot, forces_vector, dist_vector, p0=tuple(p0),
                               bounds=(bounds_min, bounds_max))
    perr = np.sqrt(np.diag(pcov))

    if linker == 'dna':
        coefficients = {'ldna': popt[0], 'pdna': popt[1], 'pprot': popt[2], 'L': popt[3:]}
        errors = {'ldna': perr[0], 'pdna': perr[1], 'pprot': perr[2], 'L': perr[3:]}
    else:
        coefficients = {'pprot': popt[0], 'L': popt[1:]}
        errors = {'pprot': perr[0], 'L': perr[1:]}

    return coefficients, errors


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


def find_close_forces(value, forces):
    return


def find_area(n, bins):
    width = bins[1]-bins[0]
    return sum(n)*width


def force_load(f, p, l, v, k=0):
    force_part = (2*p*l * (1+p*f))/(3 + 5*p*f + 8*(p*f)**5/2)
    if k != 0:
        force_part += 1/k
    return v/force_part


def fit_dudko(forces, times, v):
    force = Variable('force')
    t = Variable('t')
    gamma = Parameter('gamma', min=0.00001, value=0.1, max=1/max(forces)-0.00001)
    g = Parameter('g', min=0.00001, value=10)
    t0 = Parameter('t0', min=0.00001, value=10)
    # x = Parameter('x', min=0.1, value=1)

    # v = 1/2
    if v == 1/2:
        model = Model({t: t0 * (1-gamma*force)**(-1) * exp(-g * (1 - (1-gamma*force)**2))})
        fit = Fit(model, t=times, force=forces)
        fit_result = fit.execute()
        result = {'t0': fit_result.value(t0), 'x': 2*fit_result.value(gamma)*fit_result.value(g),
                         'g': fit_result.value(g)}

    # v = 2/3
    if v == 2/3:
        model = Model({t: t0 * (1 - gamma * force) ** (-1/2) * exp(-g * (1 - (1 - gamma * force) ** (3/2)))})
        fit = Fit(model, t=times, force=forces)
        fit_result = fit.execute()
        result = {'t0': fit_result.value(t0), 'x': (3/2) * fit_result.value(gamma) * fit_result.value(g),
                         'g': fit_result.value(g)}

    # v = 1
    if v == 1:
        model = Model({t: t0 * exp(-g * (1 - (1 - gamma * force)))})
        fit = Fit(model, t=times, force=forces)
        fit_result = fit.execute()
        result = {'t0': fit_result.value(t0), 'x': fit_result.value(gamma) * fit_result.value(g),
                       'g': fit_result.value(g)}

    return result


def fit_part_wlc(dist, forces):
    d = Variable('d')
    f = Variable('f')
    length = Parameter('length', value=max(dist) + 20, min=max(dist) + 1)
    pprot = Parameter('pprot', value=1, min=0.1)
    model = Model({f: pprot * (0.25 / ((1 - d / length) ** 2) - 0.25 + d / length)})
    fit = Fit(model, d=dist, f=forces)
    fit_result = fit.execute()
    coefficients = {'length': fit_result.value(length), 'pprot': fit_result.value(pprot)}
    return coefficients


def exp_fit(f, t0, k1, k2, v):
    return t0 * (1-k1*f)**(1-1/v) * np.exp(-k2 * (1-(1-k1*f)**v))
