import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from symfit import parameters, variables, Fit, Model, sqrt, GreaterThan, Parameter, Variable


def extract_curve(name, data_path, data_file_prefix, data_file_suffix):
    fname = data_path + data_file_prefix + name + data_file_suffix
    with open(fname, 'r') as myfile:
        # reading values
        values = []
        for line in myfile.readlines():
            if 'd' not in line:
                values.append(line.strip().split(';'))
    # extracting each trace
    for trace in range(1, len(values[0]), 2):
        pairs = [[float(values[point][trace]), float(values[point][trace + 1])]
                 for point in range(len(values)) if values[point][trace]]
        pairs = sorted(pairs)
        dist = [v[0] for v in pairs]
        force = [v[1] for v in pairs]
        yield dist, force
    return


def find_ranges(dist, gap_size, min_str_dist, break_size):
    ranges = []
    current_range = [dist[0]]
    for k in range(len(dist) - 1):
        # we are looking for the current range closure
        if len(current_range) == 1 and dist[k + 1] - dist[k] > gap_size and dist[k] - current_range[0] > min_str_dist:
            current_range.append(dist[k])
            ranges.append(current_range)
            current_range = []
        # we are looking for the current range opening
        elif len(current_range) == 0 and dist[k + 1] - dist[k] < gap_size and dist[k + 1] - ranges[-1][-1] > break_size:
            current_range = [dist[k + 1]]
    # closing the range on the final point
    if len(current_range) == 1 and dist[-1] - current_range[0] > gap_size:
        current_range.append(dist[-1])
        ranges.append(current_range)
    return ranges


def fit_curve(dist, forces, ranges, high_force_cutoff):
    dist_to_fit = [np.array([dist[_] for _ in range(len(dist))
                             if current_range[0] < dist[_] < current_range[1] and forces[_] > high_force_cutoff])
                   for current_range in ranges]
    forces_to_fit = [np.array([forces[_] for _ in range(len(dist))
                               if current_range[0] < dist[_] < current_range[1] and forces[_] > high_force_cutoff])
                     for current_range in ranges]
    # the model for symfit
    dist_variables = []
    forces_variables = []
    lenghts = []
    model_dict = {}
    elastic_moduli = Parameter('K', value=280, min=0.1)
    persistence_dna = Parameter('pd', value=0.137, min=0.01)
    persistence_protein = Parameter('pp', value=5.877, min=0.01)
    for k in range(len(ranges)):
        dist_variables.append(Variable('d' + str(k)))
        forces_variables.append(Variable('F' + str(k)))
        lenghts.append(Parameter('L' + str(k), value=max(ranges[k][1]-ranges[0][1], ranges[0][1]), min=0.1))
        if k == 0:
            model_dict[dist_variables[k]] = lenghts[0] * (1 - sqrt(persistence_dna / forces_variables[0]) + forces_variables[0] / elastic_moduli)
        else:
            model_dict[dist_variables[k]] = lenghts[0] * (1 - sqrt(persistence_dna / forces_variables[k]) + forces_variables[k] / elastic_moduli) + \
                                    lenghts[k] * (1 - sqrt(persistence_protein / forces_variables[k]))
    model = Model(model_dict)

    # fitting
    if len(ranges) == 2:
        fit = Fit(model, d0=dist_to_fit[0], d1=dist_to_fit[1], F0=forces_to_fit[0], F1=forces_to_fit[1])
    elif len(ranges) == 3:
        fit = Fit(model, d0=dist_to_fit[0], d1=dist_to_fit[1], d2=dist_to_fit[2],
                  F0=forces_to_fit[0], F1=forces_to_fit[1], F2=forces_to_fit[2])
    else:
        fit = Fit(model, d0=dist_to_fit[0], d1=dist_to_fit[1], d2=dist_to_fit[2], d3=dist_to_fit[3],
                  F0=forces_to_fit[0], F1=forces_to_fit[1], F2=forces_to_fit[2], F3=forces_to_fit[3])
    fit_result = fit.execute()

    # extracting coefficients
    coefficients = {'characteristic lengths (L)': [fit_result.value(L) for L in lenghts],
                    'DNA elastic modulus (K)': fit_result.value(elastic_moduli),
                    'pDNA/KbT': fit_result.value(persistence_dna),
                    'pProtein/KbT': fit_result.value(persistence_protein),
                    'persistence length (pProtein)': 4.114 / fit_result.value(persistence_protein),
                    'persistence length (pDNA)': 4.114 / fit_result.value(persistence_dna)}
    return coefficients
