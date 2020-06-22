import numpy as np

from symfit import Variable, Parameter, Fit, Model, exp


def merge_dicts(dict1, dict2):
    return {**dict1, **dict2}


def calc_forward_wlc(d_array, length, pprot):
    return pprot * (0.25/((1-d_array/length)**2) - 0.25 + d_array/length)


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
