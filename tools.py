import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from symfit import Fit, Model, sqrt, Parameter, Variable, exp
from scipy.stats import norm
from scipy.optimize import curve_fit

# available models for curve fitting
available_linkers = ['none', 'harmonic', 'dna']
colors = [mcolors.CSS4_COLORS['red'],
          mcolors.CSS4_COLORS['green'],
          mcolors.CSS4_COLORS['blue'],
          mcolors.CSS4_COLORS['yellow'],
          mcolors.CSS4_COLORS['cyan'],
          mcolors.CSS4_COLORS['magenta'],
          mcolors.CSS4_COLORS['orange'],
          mcolors.CSS4_COLORS['purple'],
          mcolors.CSS4_COLORS['lime']]


def extract_curves(name, data_path, data_file_prefix, data_file_suffix, unit='A'):
    dist_total = []
    forces_total = []
    fname = data_path + data_file_prefix + name + data_file_suffix
    with open(fname, 'r') as myfile:
        # reading values
        values = []
        for line in myfile.readlines():
            if 'd' not in line and 'D' not in line:
                values.append(line.strip().split(';'))
    # extracting each trace
    beg = len(values[0])%2
    end = len(values[0])
    for trace in range(beg, end, 2):
        pairs = [[float(values[point][trace]), float(values[point][trace + 1])]
                 for point in range(len(values)) if values[point][trace]]
        pairs = sorted(pairs)
        dist = [v[0] for v in pairs]
        force = [v[1] for v in pairs]
        if unit == 'nm':
            dist = [10*d for d in dist]
        dist_total.append(dist)
        forces_total.append(force)
    return dist_total, forces_total


def smooth_curves(dist, forces):
    forces_smooth = []
    errors = []
    dist_smooth = list(range(int(min(dist)), int(max(dist))+1))
    for d in dist_smooth:
        forces_interval = [forces[_] for _ in range(len(dist)) if d <= dist[_] < d + 1]
        try:
            forces_smooth.append(sum(forces_interval)/len(forces_interval))
        except ZeroDivisionError:
            errors.append(d)
    for d in list(reversed(errors)):
        dist_smooth.pop(dist_smooth.index(d))
    return dist_smooth, forces_smooth


def running_average(dist, forces):
    window = int((max(dist)-min(dist))/100)
    dist_smooth = np.convolve(dist, np.ones((window,))/window, mode='valid')
    forces_smooth = np.convolve(forces, np.ones((window,))/window, mode='valid')
    return dist_smooth, forces_smooth


def find_derivative(dist, forces):
    forces_derivative = np.diff(np.array(forces))/np.diff(np.array(dist))
    dist_derivative = np.array(dist[:-1]) + np.diff(np.array(dist))/2
    return running_average(dist_derivative, forces_derivative)


def find_ranges(dist, forces, gap=10, high_force_cutoff=5):
    dist_smooth, forces_smooth = smooth_curves(dist, forces)
    dist_smooth, forces_smooth = running_average(dist_smooth, forces_smooth)
    dist_derivative, forces_derivative = find_derivative(dist_smooth, forces_smooth)
    ranges = []
    current_range = [dist_derivative[0]]
    negative_derivative = -1
    for d in range(len(dist_derivative)):
        if forces_derivative[d] < 0:
            if negative_derivative < 0:
                negative_derivative = d
            # if the window of negative derivative is at least equal 'gap' parameter, we finish the current range
            if negative_derivative > 0 and dist_derivative[d] - dist_derivative[negative_derivative] >= gap \
                    and len(current_range) == 1 and forces_smooth[d] >= high_force_cutoff:
                current_range.append(dist_derivative[negative_derivative])
                ranges.append(current_range)
                current_range = []
        # establishing new range beginning as the first point when derivative becomes again positive
        if forces_derivative[d] > 0:
            if negative_derivative > 0 and len(current_range) == 0:
                current_range.append(dist_derivative[d])
            negative_derivative = -1
    if len(current_range) == 1:
        current_range.append(dist_derivative[-1])
        ranges.append(current_range)
    return ranges, dist_smooth, forces_smooth


def fit_last_part(dist, forces):
    dlast = Variable('dlast')
    flast = Variable('flast')
    llast = Parameter('llast', value=max(dist) + 20, min=max(dist) + 1)
    plast = Parameter('plast', value=1, min=0.1)
    model = Model({flast: plast * (0.25 / ((1 - dlast / llast) ** 2) - 0.25 + dlast / llast)})
    fit = Fit(model, dlast=dist, flast=forces)
    fit_result = fit.execute()
    coefficients = {'llast': fit_result.value(llast), 'plast': fit_result.value(plast)}
    return coefficients


def solve_wlc(f, p, l=None, x=None):
    b = -2.25 - f/p
    c = 1.5 + 2*f/p
    d = - f/p
    roots = np.roots([1, b, c, d])
    result = 0
    for root in roots[np.imag(roots) == 0]:
        if 0 < root < 1:
            result = np.real(root)
            break
    if l:
        return l*result
    elif x:
        return x/result
    else:
        return 0


def separate_linker(ranges, dist_smooth, forces_smooth, linker='none'):
    coefficients = {}
    if linker == 'none':
        return ranges, dist_smooth, forces_smooth, coefficients
    elif linker == 'harmonic':
        # plt.plot(dist_smooth,forces_smooth)
        dist = [dist_smooth[_] for _ in range(len(dist_smooth)) if ranges[-1][0] <= dist_smooth[_] <= ranges[-1][1]]
        forces = [forces_smooth[_] for _ in range(len(dist_smooth)) if ranges[-1][0] <= dist_smooth[_] <= ranges[-1][1]]
        dlinker = Variable('dlinker')
        flinker = Variable('flinker')
        plinker = Parameter('plinker', min=0.1, value=0.3)
        klinker = Parameter('klinker', min=0, value=1)
        plt.plot(dist, forces, label='original', linewidth=2)
        for k in np.linspace(1, 100, 10):
            print(k)
            dist_protein = [d-f/k for d,f in zip(dist, forces) if f > 6]
            forces_to_plot = [f for f in forces if f > 6]
            llinker = Parameter('llinker', value=max(dist_protein) + 20, min=max(dist_protein) + 1)
            model = Model({flinker: plinker * (0.25 / ((1 - dlinker / llinker) ** 2) - 0.25 + dlinker / llinker)})
            fit = Fit(model, dlinker=dist_protein, flinker=forces_to_plot)
            fit_result = fit.execute()
            coefficients['llinker'] = (fit_result.value(llinker))
            coefficients['plinker'] = (fit_result.value(plinker))
            print(coefficients)
            forces_protein = coefficients['plinker'] * (0.25 / ((1 - np.array(dist_protein) / coefficients['llinker']) ** 2) - 0.25 + np.array(dist_protein) / coefficients['llinker'])
            # plt.plot(dist_protein, forces_protein, label=str(k))
            dist_all = dist_protein + forces_protein/k
            plt.plot(dist_all, forces_protein, label=str(k))
        plt.legend()
        plt.show()
        return ranges, dist_smooth, forces_smooth, coefficients
    elif linker == 'dna':
        low_force_cutoff = 0.1
        dist_to_fit = [dist_smooth[_] for _ in range(len(dist_smooth)) if
                       ranges[0][0] <= dist_smooth[_] <= ranges[0][1] and forces_smooth[_] > low_force_cutoff]
        forces_to_fit = [forces_smooth[_] for _ in range(len(dist_smooth)) if
                       ranges[0][0] <= dist_smooth[_] <= ranges[0][1] and forces_smooth[_] > low_force_cutoff]
        ddna = Variable('ddna')
        fdna = Variable('fdna')
        ldna = Parameter('ldna', value=max(dist_to_fit) + 20, min=max(dist_to_fit) + 1)
        pdna = Parameter('pdna', value=1, min=0.1)
        model = Model({fdna: pdna * (0.25 / ((1 - ddna / ldna) ** 2) - 0.25 + ddna / ldna)})
        fit = Fit(model, ddna=dist_to_fit, fdna=forces_to_fit)
        fit_result = fit.execute()
        coefficients['ldna'] = (fit_result.value(ldna))
        coefficients['pdna'] = (fit_result.value(pdna))
        dist_dna = [solve_wlc(f, coefficients['pdna'], l=coefficients['ldna']) for f in forces_smooth]
        dist_protein = [d0-d1 for d0, d1 in zip(dist_smooth, dist_dna) if d0 > ranges[1][0]]
        forces_protein = [forces_smooth[_] for _ in range(len(forces_smooth)) if dist_smooth[_] > ranges[1][0]]
        return ranges, dist_protein, forces_protein, coefficients
    else:
        raise NameError("Unknown linker for analysis. Known models are " + str(available_linkers) + '.')


def fit_curve_dna(ranges, dist, forces):
    low_force_cutoff = 0.1
    dist_to_fit = [dist[_] for _ in range(len(dist)) if
                   ranges[0][0] <= dist[_] <= ranges[0][1] and forces[_] > low_force_cutoff]
    forces_to_fit = [forces[_] for _ in range(len(dist)) if
                     ranges[0][0] <= dist[_] <= ranges[0][1] and forces[_] > low_force_cutoff]
    ddna = Variable('ddna')
    fdna = Variable('fdna')
    ldna = Parameter('ldna', value=max(dist_to_fit) + 20, min=max(dist_to_fit) + 1)
    pdna = Parameter('pdna', value=1, min=0.1)
    model = Model({fdna: pdna * (0.25 / ((1 - ddna / ldna) ** 2) - 0.25 + ddna / ldna)})
    fit = Fit(model, ddna=dist_to_fit, fdna=forces_to_fit)
    fit_result = fit.execute()
    coefficients = {'ldna': fit_result.value(ldna), 'pdna': fit_result.value(pdna)}

    # dist_dna = [np.array([solve_wlc(f, coefficients['pdna'], l=coefficients['ldna']) for f in force_range])
    #                 for force_range in forces_to_fit]
    dist_dna_smooth = np.array([solve_wlc(f, coefficients['pdna'], l=coefficients['ldna']) for f in forces])
    # plt.plot(forces, dist_dna_smooth)
    # plt.plot(forces, dist)
    dist_protein = [d0-d1 for d0, d1 in zip(dist, dist_dna_smooth) if d0 >= ranges[1][0]]
    forces_protein = [forces[_] for _ in range(len(forces)) if dist[_] >= ranges[1][0]]
    new_ranges = [[0,200], [400,700]]
    dist_protein_to_fit = [np.array([dist_protein[_] for _ in range(len(dist_protein))
                             if current_range[0] <= dist_protein[_] <= current_range[1] and
                             forces_protein[_] > low_force_cutoff]) for current_range in new_ranges]
    forces_protein_to_fit = [np.array([forces_protein[_] for _ in range(len(dist_protein))
                             if current_range[0] <= dist_protein[_] <= current_range[1] and
                               forces_protein[_] > low_force_cutoff]) for current_range in new_ranges]
    coefficients = fit_last_part(dist_protein_to_fit[-1], forces_protein_to_fit[-1])
    print(coefficients)
    plt.plot(dist_protein, forces_protein)
    plt.show()
    return coefficients


def fit_curve(ranges, dist, forces, temperature=None, elastic_moduli_boudaries=(250, 400)):
    low_force_cutoff = 0.1
    if temperature:
        temperature_factor = temperature * 0.0138
    else:
        temperature_factor = 4.114

    dist_to_fit = [np.array([dist[_] for _ in range(len(dist))
                             if current_range[0] <= dist[_] <= current_range[1] and
                             forces[_] > low_force_cutoff]) for current_range in ranges]
    forces_to_fit = [np.array([forces[_] for _ in range(len(dist))
                             if current_range[0] <= dist[_] <= current_range[1] and
                               forces[_] > low_force_cutoff]) for current_range in ranges]

    coefficients = fit_last_part(dist_to_fit[-1], forces_to_fit[-1])

    # the model for symfit
    dist_variables = []
    forces_variables = []
    lengths = []
    persistence_protein = Parameter('pp', value=0.3, min=0)  # , min=3, max=10)
    model_dict = {}

    for k in range(len(ranges)):
        dist_variables.append(Variable('d' + str(k)))
        forces_variables.append(Variable('F' + str(k)))
        lengths.append(Parameter('Lp' + str(k), value=max(dist_to_fit[k])+20, min=max(dist_to_fit[k])+1,
                       max=coefficients['llast']))
        model_dict[forces_variables[k]] = persistence_protein * (
                    0.25 / ((1 - dist_variables[k] / lengths[k]) ** 2) - 0.25 + dist_variables[k]/lengths[k])

    model = Model(model_dict)
    arguments = {}
    for k in range(len(ranges)):
        arguments['d' + str(k)] = dist_to_fit[k]
        arguments['F' + str(k)] = forces_to_fit[k]

    fit = Fit(model, **arguments)
    fit_result = fit.execute()

    # extracting coefficients
    coefficients['L'] = [fit_result.value(L) for L in lengths]
    coefficients['pProtein/KbT'] = fit_result.value(persistence_protein)
    coefficients['pProtein'] = temperature_factor / fit_result.value(persistence_protein)
    return coefficients


def find_rupture_forces(dist, forces, ranges, coefficients):
    lengths = coefficients['L']
    rupture_forces = {}
    for r in range(len(ranges)):
        current_range = ranges[r]
        length = lengths[r]
        rupture_force = max([forces[_] for _ in range(len(dist)) if current_range[0] < dist[_] < current_range[1]])
        rupture_forces[length] = rupture_force
    return rupture_forces


def transform_coordinates(dist, forces, ranges, coefficients, low_force_cutoff=0.1):
    result = []
    p = coefficients['pProtein/KbT']
    for r in range(len(ranges)):
        current_range = ranges[r]
        dist_part = [dist[_] for _ in range(len(dist)) if current_range[0] <= dist[_] <= current_range[1]
                     and forces[_] > low_force_cutoff]
        forces_part = [forces[_] for _ in range(len(dist)) if current_range[0] <= dist[_] <= current_range[1]
                     and forces[_] > low_force_cutoff]
        result += [solve_wlc(forces_part[_], p, x=dist_part[_]) for _ in range(len(dist_part))]
    return result


def calculate_protein_force(d, L, p):
    return p * (0.25/((1-d/L)**2) - 0.25 + d/L)


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


def make_histograms(coefficients_total, contour_lengths, rupture_forces, name, residues_distance=3.88, cluster_max_gap=15,
                    show_plots=False):
    fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))
    results = {}

    # contour length histogram
    axes[0, 0].set_title('Contour length histogram')
    axes[0, 0].set_xlabel('Countour length [A]')
    axes[0, 0].set_ylabel('Probability')
    all_coefficients = [coefficient for r in coefficients_total for coefficient in r['L']]
    coefficients_clusters = cluster_coefficients(all_coefficients, maxgap=cluster_max_gap)
    min_l = min(all_coefficients)
    max_l = max(all_coefficients)
    axes[0, 0].set_xlim(min_l, max_l)
    x = np.linspace(min_l, max_l, 100)
    results['contour_length_histo'] = []
    for r in range(len(coefficients_clusters)-1):
        n, bins, patches = axes[0, 0].hist(coefficients_clusters[r], density=True, facecolor=colors[r % len(colors)],
                                           alpha=0.5)
        (mu, sigma) = norm.fit(coefficients_clusters[r])
        residues = int(round(mu / residues_distance, 0)) + 1
        results['contour_length_histo'].append([round(mu, 3), round(sigma, 3), residues])
        label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
        axes[0, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
    n, bins, patches = axes[0, 0].hist(coefficients_clusters[-1], bins=np.arange(min_l, max_l + 5, 5),
                                           density=True, facecolor='k', alpha=0.5)
    axes[0, 0].legend()

    # forces histogram
    axes[1, 0].set_title('Rupture forces histogram')
    axes[1, 0].set_xlabel('Rupture force [pN]')
    axes[1, 0].set_ylabel('Probability')
    min_f = min(list(rupture_forces.values()))
    max_f = max(list(rupture_forces.values()))
    axes[1, 0].set_xlim(min_f, max_f)
    x = np.linspace(min_f, max_f, 100)
    results['rupture_forces'] = []
    dudko_analysis = []
    for r in range(len(coefficients_clusters)-1):
        force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[r]]
        n, bins, patches = axes[1, 0].hist(force_cluster, density=True, facecolor=colors[r % len(colors)], alpha=0.5)
        dudko_analysis.append([n, bins])
        (mu, sigma) = norm.fit(force_cluster)
        results['rupture_forces'].append([round(mu, 3), round(sigma, 3)])
        label = str(round(mu, 3)) + ' (' + str(results['contour_length_histo'][r][0]) + ')'
        axes[1, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
    force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[-1]]
    n, bins, patches = axes[1, 0].hist(force_cluster, bins=np.arange(min_l, max_l + 5, 5),
                                           density=True, facecolor='k', alpha=0.5)
    axes[1, 0].legend()

    # transformed contour length histogram
    axes[0, 1].set_title('Transformed contour length histogram')
    axes[0, 1].set_xlabel('Contour length [pN]')
    axes[0, 1].set_ylabel('Probability')
    length_clusters = cluster_coefficients(contour_lengths, maxgap=1, minspan=10)
    min_l = min([length for cluster in length_clusters[:-1] for length in cluster])
    max_l = max([length for cluster in length_clusters[:-1] for length in cluster])
    axes[0, 1].set_xlim(min_l, max_l)
    x = np.linspace(min_l, max_l, 100)
    results['transformed'] = []
    for r in range(len(length_clusters)-1):
        n, bins, patches = axes[0, 1].hist(length_clusters[r], density=True,
                                           facecolor=colors[r % len(colors)], alpha=0.5)
        (mu, sigma) = norm.fit(length_clusters[r])
        residues = int(round(mu / residues_distance, 0)) + 1
        results['transformed'].append([round(mu, 3), round(sigma, 3), residues])
        label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
        axes[0, 1].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
    if len(length_clusters[-1]) > 4:
        n, bins, patches = axes[0, 0].hist(length_clusters[-1], bins=np.arange(min_l, max_l + 5, 5),
                                           density=True, facecolor='k', alpha=0.5)
    axes[0, 1].legend()

    # Dudko analysis
    axes[1, 1].set_title('Dudko analysis')
    axes[1, 1].set_xlabel('Rupture force [pN]')
    axes[1, 1].set_ylabel('Loading rate')

    # axes[1, 1].legend()


    fig.tight_layout()
    if show_plots:
        plt.show()
    else:
        print("Saving histograms figure to " + name + '_histograms.png')
        plt.savefig(name + '_histograms.png')
    return results


def make_partial_plots(dist_total, forces_total, ranges_total, coefficients_total, name, linker='None', show_plots=False,
                      columns=4, residues_distance=3.88):
    rows = max(int(np.ceil(float(len(ranges_total))/columns)), 2)
    fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
    results = []
    for k in range(len(ranges_total)):
        row = int(k/4)
        col = k % 4
        dist = dist_total[k]
        forces = forces_total[k]
        ranges = ranges_total[k]
        coefficients = coefficients_total[k]
        axes[row, col].set_xlim(min(dist), max(dist))
        axes[row, col].set_ylim(0, max(forces))
        axes[row, col].set_title(name + '_' + str(k+1))
        axes[row, col].set_xlabel('Extension [A]')
        axes[row, col].set_ylabel('Force [pN]')
        axes[row, col].plot(np.array(dist), np.array(forces))
        for l in range(len(coefficients['L'])):
            length = coefficients['L'][l]
            residues = int(round(length / residues_distance, 0)) + 1
            label = str(round(length, 3)) + ' (' + str(residues) + ' AA)'
            d_plot = np.linspace(min(dist), coefficients['L'][l]-1)
            forces_plot = calculate_protein_force(d_plot, coefficients['L'][l], coefficients['pProtein/KbT'])
            axes[row, col].plot(d_plot, forces_plot, label=label)
        for r in ranges:
            axes[row, col].axvline(x=r[0], linestyle='--', color='#808080', linewidth=0.5)
            axes[row, col].axvline(x=r[1], linestyle='--', color='#808080', linewidth=0.5)
        axes[row, col].legend()
    fig.tight_layout()
    if show_plots:
        plt.show()
    else:
        print("Saving contour lengths figure to " + name + '_contour_lengths.png')
        plt.savefig(name + '_contour_lengths.png')
    return results


def find_energies(dist, forces, coefficients):
    result = 0
    return result


def exp_fit(f, t0, k1, k2, v):
    return t0 * (1-k1*f)**(1-1/v) * np.exp(-k2 * (1-(1-k1*f)**v))


def find_averages(dist_total, forces_total):
    dist_aver = []
    forces_aver = []
    return


def find_force_load(f, aver_dist, aver_forces, extension_speed):
    return 1


def extract_life_times(force_counts, force_bins, aver_dist, aver_forecs, extension_speed, show_plot):
    width = force_bins[1]-force_bins[0]

    curve = []
    forces = []
    for k in range(1, len(force_bins)):
        forces.append(force_bins[0] + (float(k)/2) * width)
        curve.append(((force_counts[k-1]/2 + sum(force_counts[k:]))*width) /
                     (force_counts[k-1] * find_force_load(forces[-1], aver_dist, aver_forecs, extension_speed)))
    popt, pcov = curve_fit(exp_fit, forces, curve, p0=(1, 0.01, 0.1, 0.5))
    print(popt)
    # x = np.linspace(min(forces), max(forces), 100)
    # y = exp_fit(forces, popt[0], popt[1], popt[2], popt[3])
    plt.plot(forces, np.array(curve), 'ro')
    if show_plot:
        plt.show()
        plt.close('all')
    else:
        plt.savefig('Dudko_plot.png')
    xd = 0
    g = 0
    v = 0
    tau = 0

    return xd, g, v, tau


def save_data(name, details, ranges, coefficients, rupture_forces, contour_length_gain, histograms):
    fname = name + '_results'
    result = []
    # intro
    result.append('Analyzed file: ' + name)
    result.append('Number of residues: ' + str(details['residues']))
    result.append('Distance between the termini: ' + str(details['distance']) + 'A')
    result.append('Linker to tweezers: ' + str(details['linker']))
    result.append('Unit of distance data: ' + str(details['unit']))
    result.append('Source of the data: ' + details['source'])
    result.append('----')

    # parameters
    result.append('Parameters used for the calculation:')
    result.append('----')

    # ranges
    result.append('Found ' + str(len(ranges)) + ' traces.')
    result.append('The ranges between the jumps are:')
    for k in range(len(ranges)):
        result.append(str(k+1) + '\t' + ';'.join([str(round(x, 3)) + '-' + str(round(y, 3)) for x, y in ranges[k]]))
    result.append('----')

    # coefficients
    result.append('Coefficients obtained:')
    result.append('trace\tcoefficient\tvalue')
    for k in range(len(coefficients)):
        for key in coefficients[k].keys():
            result.append(str(k) + '\t' + key + ':\t' + str(coefficients[k][key]))
    result.append('----')

    # rupture forces
    result.append('Rupture forces:')
    result.append('contour length\tvalue')
    for key in rupture_forces.keys():
        result.append(str(key) + '\t' + str(rupture_forces[key]))
    result.append('----')

    # histogram analysis
    result.append('The contour lengths found by fitting histograms:')
    for length in histograms['contour_length_histo']:
        result.append('length: ' + str(length[0]) + '\tsigma: ' + str(length[1]) + '\t' + str(length[2]) + ' residues')
    result.append('----')
    result.append('Rupture forces histogram:')
    for k in range(len(histograms['rupture_forces'])):
        value = histograms['rupture_forces'][k]
        length = histograms['contour_length_histo'][k][0]
        result.append('force: ' + str(value[0]) + '\tsigma: ' + str(value[1]) + '\tlength: ' + str(length))
    result.append('----')


    # saving to file
    with open(fname, 'w') as myfile:
        myfile.write('\n'.join(result))
    print('The results saved to ' + fname)
    return


def merge_dicts(dict1, dict2):
    return {**dict1, **dict2}