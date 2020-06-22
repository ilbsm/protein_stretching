import numpy as np
from matplotlib import pyplot as plt
from stretchme.tools import invert_wlc_np
from scipy.stats import cauchy, ks_2samp

parameters = {
    'force_range': [1, 40],                 # the range of stretching forces
    'number': 20,                           # number of spectra to simulate
    'force_variance': 1,                  # the variance of blur of applied force
    'dist_variance': 1,                   # the variance of blur of measured extension
    'persistence': 5.88,                      # protein persistence length/kbT
    'length': 80,                          # protein contour length
    'length_dna': 360,                     # DNA linker contour length
    'persistence_dna': 0.16,                 # DNA persistence length/kbT
    'add_dna': True,                        # if to add the DNA linker
    'low_force_cutoff': 1                 # the cutoff of reliable force detection
}


def calc_dist(f_array, length, p):
    return length * np.array([invert_wlc_np(f, p) for f in f_array])


def simulate_spectrum():
    f_space = np.linspace(parameters['force_range'][0], parameters['force_range'][1], 1000)
    forces_mixed = f_space + np.random.normal(0, parameters['force_variance'], 1000)
    dist_exact = calc_dist(f_space, parameters['length'], parameters['persistence'])
    dist_mixed = dist_exact + np.random.normal(0, parameters['dist_variance'], 1000)
    if parameters['add_dna']:
        dist_exact_dna = calc_dist(f_space, parameters['length_dna'], parameters['persistence_dna'])
        dist_mixed_dna = dist_exact_dna + np.random.normal(0, parameters['dist_variance'], 1000)
        dist_mixed += dist_mixed_dna
    trace = [[d, f] for d, f in zip(dist_mixed, forces_mixed)]
    trace.sort(key=lambda x: x[0])
    return trace


def simulate_experiment():
    traces = [simulate_spectrum() for k in range(parameters['number'])]
    return traces


def calcualate_contour_lengths(data, pprot):
    contour_lengths = []
    for trace in data:
        dist = np.array([d for d, f in trace if f > parameters['low_force_cutoff']])
        xs = np.array([invert_wlc_np(f, pprot) for d, f in trace if
                       f > parameters['low_force_cutoff']])
        contour_lengths += list(dist / xs)
    return contour_lengths


def calculate_pvalue(data, pprot):
    contour_lengths = calcualate_contour_lengths(data, pprot)
    mu, gamma = cauchy.fit(contour_lengths)
    test_statistic = cauchy.rvs(mu, gamma, len(contour_lengths))
    result = [mu, gamma, pprot, ks_2samp(contour_lengths, test_statistic).pvalue]
    return result


def plot_results(data, parameters):
    mu, gamma, pprot, pvalue = parameters
    contour_lengths = calcualate_contour_lengths(data, pprot)
    lengths = np.linspace(min(contour_lengths), max(contour_lengths))
    label = 'Cauchy(' + str(round(mu, 3)) + ', ' + str(round(gamma, 3)) + ')'
    plt.title('pProt=' + str(round(pprot, 3)) + '; Pvalue=' + str(round(pvalue, 3)))
    plt.hist(contour_lengths, bins=200, density=True)
    plt.plot(lengths, cauchy.pdf(lengths, mu, gamma), linestyle='--', label=label)
    plt.xlim(min(contour_lengths), max(contour_lengths))
    plt.legend()
    plt.show()
    return


def find_parameters(data, precision=0.0001, initial=0.3):
    pmin = 0
    pmax = 0
    results = {}
    # looking for initial range
    found = False
    for pprot in np.linspace(initial, 50.3, 101):
        res = calculate_pvalue(data, pprot)[-1]
        if res > 0.00001:
            found = True
            pmin = pprot - (50.3 - initial)/100
            results[pmin] = calculate_pvalue(data, pmin)[-1]
        if found and res < 0.00001:
            pmax = pprot
            results[pmax] = res
            break

    # looking in the interval - the golden division algorithm
    p0 = pmin
    p3 = pmax
    epsilon = p3 - p0
    while epsilon > precision:
        p1 = 2*p0/3 + p3/3
        p2 = p0/3 + 2*p3/3
        results[p1] = calculate_pvalue(data, p1)[-1]
        results[p2] = calculate_pvalue(data, p2)[-1]
        if results[p1] > results[p2]:
            p3 = p2
        else:
            p0 = p1
        epsilon = p3 - p0
    pprot = (p3+p0)/2
    result = calculate_pvalue(data, pprot)
    plot_results(data, result)
    return result


data = simulate_experiment()
x1 = []
x2 = []
res = []
for trace in data:
    x1 += [invert_wlc_np(f, parameters['persistence']) for d, f in trace if f > parameters['low_force_cutoff']]
    x2 += [invert_wlc_np(f, parameters['persistence_dna']) for d, f in trace if f > parameters['low_force_cutoff']]
    res += [d for d, f in trace if f > parameters['low_force_cutoff']]
values = res/(np.array(x1)+4.5*np.array(x2))
n, bins, patches = plt.hist(values, bins=200, density=True)
f_space = np.linspace(min(bins), max(bins))
mu, gamma = cauchy.fit(values)
plt.plot(f_space, cauchy.pdf(f_space, mu, gamma), label=str(round(mu, 3)))
test_statistic = cauchy.rvs(mu, gamma, len(res))
pvalue = ks_2samp(values, test_statistic).pvalue
plt.legend()
plt.show()
print(pvalue)
# result = find_parameters(data)


''' Main part '''
