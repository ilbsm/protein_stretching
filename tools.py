import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from symfit import Variable, Parameter, Fit, Model, exp


def merge_dicts(dict1, dict2):
    return {**dict1, **dict2}


def invert_wlc_np(f, p):
    if f == 0:
        return 0
    coefs = [1, -(2.25 + f / p), (1.5 + 2 * f / p), -f / p]
    result = np.roots(coefs)
    result = np.real(result[np.isreal(result)])
    result = result[result < 1]
    result = result[result > 0]
    if len(result) == 1:
        return result[0]
    else:
        return 0


def running_average(x, y):
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
    clusters.append([loc for k in singles for loc in clusters[k]])
    for k in list(reversed(singles)):
        clusters.pop(k)
    return clusters


def find_hist_ranges(data):
    n, bins, patches = plt.hist(data, bins=200, density=True, alpha=0.0)
    xs = []
    for k in range(len(bins) - 1):
        xs.append((bins[k] + bins[k + 1]) / 2)

    diff_x = np.diff(xs)
    diff_n = np.diff(n)

    to_cluster = []
    for x, y in zip(xs[:-1], diff_n / diff_x):
        if abs(y) < 0.00001:
            to_cluster.append(x)

    clusters = cluster_coefficients(to_cluster, minnumber=2)
    fit_ranges = [[min(bins)]]
    for k in range(len(clusters) - 1):
        cluster = clusters[k]
        pos = (max(cluster) + min(cluster)) / 2
        fit_ranges[-1].append(pos)
        fit_ranges.append([pos])
    fit_ranges.pop(-1)
    return fit_ranges


def force_load(f, p, l, v, k=0):
    force_part = (2*p*l * (1+p*f))/(3 + 5*p*f + 8*(p*f)**5/2)
    if k != 0:
        force_part += 1/k
    return v/force_part


def fit_dudko(forces, times):
    force = Variable('force')
    t = Variable('t')
    gamma = Parameter('gamma', min=0.00001, value=0.1, max=1/max(forces)-0.00001)
    g = Parameter('g', min=0.00001, value=10)
    t0 = Parameter('t0', min=0.00001, value=10)
    # x = Parameter('x', min=0.1, value=1)
    result = {}

    # v = 1/2
    model = Model({t: t0 * (1-gamma*force)**(-1) * exp(-g * (1 - (1-gamma*force)**2))})
    fit = Fit(model, t=times, force=forces)
    fit_result = fit.execute()
    result['1/2'] = {'t0': fit_result.value(t0), 'x': 2*fit_result.value(gamma)*fit_result.value(g),
                     'g': fit_result.value(g)}

    # v = 2/3
    model = Model({t: t0 * (1 - gamma * force) ** (-1/2) * exp(-g * (1 - (1 - gamma * force) ** (3/2)))})
    fit = Fit(model, t=times, force=forces)
    fit_result = fit.execute()
    result['2/3'] = {'t0': fit_result.value(t0), 'x': (3/2) * fit_result.value(gamma) * fit_result.value(g),
                     'g': fit_result.value(g)}

    # v = 1
    model = Model({t: t0 * exp(-g * (1 - (1 - gamma * force)))})
    fit = Fit(model, t=times, force=forces)
    fit_result = fit.execute()
    result['1'] = {'t0': fit_result.value(t0), 'x': fit_result.value(gamma) * fit_result.value(g),
                   'g': fit_result.value(g)}

    return result

##############

#
#
# def fit_last_part(dist, forces):
#     dlast = Variable('dlast')
#     flast = Variable('flast')
#     llast = Parameter('llast', value=max(dist) + 20, min=max(dist) + 1)
#     plast = Parameter('plast', value=1, min=0.1)
#     model = Model({flast: plast * (0.25 / ((1 - dlast / llast) ** 2) - 0.25 + dlast / llast)})
#     fit = Fit(model, dlast=dist, flast=forces)
#     fit_result = fit.execute()
#     coefficients = {'llast': fit_result.value(llast), 'plast': fit_result.value(plast)}
#     return coefficients
#
#
# def solve_wlc(f, p, l=None, x=None):
#     b = -2.25 - f/p
#     c = 1.5 + 2*f/p
#     d = - f/p
#     roots = np.roots([1, b, c, d])
#     result = 0
#     for root in roots[np.imag(roots) == 0]:
#         if 0 < root < 1:
#             result = np.real(root)
#             break
#     if l:
#         return l*result
#     elif x:
#         return x/result
#     else:
#         return 0
#
#
# def fit_curve(ranges, dist, forces, linker, temperature=None, elastic_moduli_boudaries=(250, 400)):
#     if linker == 'dna':
#         return fit_curve_dna(ranges, dist, forces, temperature=temperature)
#     elif linker == 'none':
#         return fit_curve_protein(ranges, dist, forces, temperature=temperature)
#     else:
#         raise NameError("Unknown linker for analysis. Known models are " + str(available_linkers) + '.')
#
#
# def fit_curve_dna(ranges, dist, forces, temperature=None):
#     low_force_cutoff = 0.1
#     if temperature:
#         temperature_factor = temperature * 0.0138
#     else:
#         temperature_factor = 4.114
#
#     dist_to_fit = [np.array([dist[_] for _ in range(len(dist))
#                              if current_range[0] <= dist[_] <= current_range[1] and
#                              forces[_] > low_force_cutoff]) for current_range in ranges]
#     forces_to_fit = [np.array([forces[_] for _ in range(len(dist))
#                              if current_range[0] <= dist[_] <= current_range[1] and
#                                forces[_] > low_force_cutoff]) for current_range in ranges]
#     for r in dist_to_fit:
#         print(str(list(r)))
#     print(';;;')
#     for r in forces_to_fit:
#         print(str(list(r)))
#
#     coefficients = {}
#
#     # the model for symfit
#     ds = []
#     fs = []
#     lengths = []
#     p_prot = Parameter('pp', value=0.3, min=0.001)  # , min=3, max=10)
#     p_dna = Parameter('pdna', value=0.3, min=0.001)  # , min=3, max=10)
#     ldna = Parameter('Ldna', value=max(dist_to_fit[0])+20, min=max(dist_to_fit[0])+1, max=5000)
#     model_dict = {}
#
#     for k in range(len(ranges)):
#         ds.append(Variable('d' + str(k)))
#         fs.append(Variable('F' + str(k)))
#         lengths.append(Parameter('Lp' + str(k), value=300, min=10, max=1000))
#         model_dict[ds[k]] = lengths[k] * invert_wlc(fs[k], p_prot)
#                             # + \
#                             # 2*ldna * solveset(x**3 - (2.25+fs[k]/p_dna)*x**2 + (1.5+2*fs[k]/p_dna)*x -
#                             #                   fs[k]/p_dna, x, domain=Interval(0, 1))
#
#     model = Model(model_dict)
#     arguments = {}
#     for k in range(len(ranges)):
#         arguments['d' + str(k)] = dist_to_fit[k]
#         arguments['F' + str(k)] = forces_to_fit[k]
#
#     fit = Fit(model, **arguments)
#     fit_result = fit.execute()
#
#     # extracting coefficients
#     coefficients['LDNA'] = fit_result.value(ldna)
#     coefficients['pDNA/KbT'] = fit_result.value(p_dna)
#     coefficients['L'] = [fit_result.value(L) for L in lengths]
#     coefficients['pProtein/KbT'] = fit_result.value(p_prot)
#     coefficients['pProtein'] = temperature_factor / fit_result.value(p_prot)
#     print(coefficients)
#     return coefficients
#
#
# def fit_curve_protein(ranges, dist, forces, temperature=None):
#     low_force_cutoff = 0.1
#     if temperature:
#         temperature_factor = temperature * 0.0138
#     else:
#         temperature_factor = 4.114
#
#     dist_to_fit = [np.array([dist[_] for _ in range(len(dist))
#                              if current_range[0] <= dist[_] <= current_range[1] and
#                              forces[_] > low_force_cutoff]) for current_range in ranges]
#     forces_to_fit = [np.array([forces[_] for _ in range(len(dist))
#                              if current_range[0] <= dist[_] <= current_range[1] and
#                                forces[_] > low_force_cutoff]) for current_range in ranges]
#
#     coefficients = fit_last_part(dist_to_fit[-1], forces_to_fit[-1])
#
#     # the model for symfit
#     dist_variables = []
#     forces_variables = []
#     lengths = []
#     persistence_protein = Parameter('pp', value=0.3, min=0.001)  # , min=3, max=10)
#     model_dict = {}
#
#     for k in range(len(ranges)):
#         dist_variables.append(Variable('d' + str(k)))
#         forces_variables.append(Variable('F' + str(k)))
#         lengths.append(Parameter('Lp' + str(k), value=max(dist_to_fit[k])+20, min=max(dist_to_fit[k])+1,
#                        max=coefficients['llast']))
#         model_dict[forces_variables[k]] = persistence_protein * (
#                     0.25 / ((1 - dist_variables[k] / lengths[k]) ** 2) - 0.25 + dist_variables[k]/lengths[k])
#
#     model = Model(model_dict)
#     arguments = {}
#     for k in range(len(ranges)):
#         arguments['d' + str(k)] = dist_to_fit[k]
#         arguments['F' + str(k)] = forces_to_fit[k]
#
#     fit = Fit(model, **arguments)
#     fit_result = fit.execute()
#
#     # extracting coefficients
#     coefficients['L'] = [fit_result.value(L) for L in lengths]
#     coefficients['pProtein/KbT'] = fit_result.value(persistence_protein)
#     coefficients['pProtein'] = temperature_factor / fit_result.value(persistence_protein)
#     return coefficients
#
#
#
#
# def transform_coordinates(dist, forces, ranges, coefficients, linker, low_force_cutoff=0.1):
#     if linker == 'none':
#         return transform_coordinates_protein(dist, forces, ranges, coefficients, low_force_cutoff=low_force_cutoff)
#     elif linker == 'dna':
#         return transform_coordinates_dna(dist, forces, ranges, coefficients, low_force_cutoff=low_force_cutoff)
#     else:
#         raise NameError("Unknown linker for analysis. Known models are " + str(available_linkers) + '.')
#
#
# def transform_coordinates_protein(dist, forces, ranges, coefficients, low_force_cutoff=0.1):
#     result = []
#     p = coefficients['pProtein/KbT']
#     for r in range(len(ranges)):
#         current_range = ranges[r]
#         dist_part = [dist[_] for _ in range(len(dist)) if current_range[0] <= dist[_] <= current_range[1]
#                      and forces[_] > low_force_cutoff]
#         forces_part = [forces[_] for _ in range(len(dist)) if current_range[0] <= dist[_] <= current_range[1]
#                      and forces[_] > low_force_cutoff]
#         result += [solve_wlc(forces_part[_], p, x=dist_part[_]) for _ in range(len(dist_part))]
#     return result
#
#
# def transform_coordinates_dna(dist, forces, ranges, coefficients, low_force_cutoff=0.1):
#     result = []
#     p = coefficients['pProtein/KbT']
#     pdna = coefficients['pDNA/KbT']
#     (mu, sigma) = norm.fit([solve_wlc(forces[_], p, x=dist[_]) for _ in range(len(dist)) if
#                             ranges[0][0] <= dist[_] <= ranges[0][1] and forces[_] > low_force_cutoff])
#     ldna = mu
#     for r in range(1, len(ranges)):
#         current_range = ranges[r]
#         forces_part = [forces[_] for _ in range(len(dist)) if current_range[0] <= dist[_] <= current_range[1]
#                        and forces[_] > low_force_cutoff]
#         dist_part = [dist[_] - solve_wlc(forces[_], pdna, l=ldna) for _ in range(len(dist))
#                      if current_range[0] <= dist[_] <= current_range[1] and forces[_] > low_force_cutoff]
#         result += [solve_wlc(forces_part[_], p, x=dist_part[_]) for _ in range(len(dist_part))]
#     return result
#
#
# def calculate_protein_force(d, L, p):
#     return p * (0.25/((1-d/L)**2) - 0.25 + d/L)
#
#
# def cluster_coefficients(coefficients, maxgap=15, minnumber=4, minspan=-1):
#     coefficients.sort()
#     clusters = [[coefficients[0]]]
#     for x in coefficients[1:]:
#         if abs(x - clusters[-1][-1]) <= maxgap:
#             clusters[-1].append(x)
#         else:
#             clusters.append([x])
#     singles = [k for k in range(len(clusters)) if len(clusters[k]) <= minnumber or
#                (max(clusters[k])-min(clusters[k])) < minspan]
#     clusters.append([loc for k in singles for loc in clusters[k]])
#     for k in list(reversed(singles)):
#         clusters.pop(k)
#     return clusters
#
#
#
#
#
# def exp_fit(f, t0, k1, k2, v):
#     return t0 * (1-k1*f)**(1-1/v) * np.exp(-k2 * (1-(1-k1*f)**v))
#
#
