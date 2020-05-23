import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from symfit import Fit, Model, sqrt, Parameter, Variable, exp
from scipy.stats import norm
from scipy.optimize import curve_fit

# available models for curve fitting
available_linkers = ['none', 'harmonic', 'dna']


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
        if unit == 'np':
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
    window =int((max(dist)-min(dist))/100)
    dist_smooth = np.convolve(dist, np.ones((window,))/window, mode='valid')
    forces_smooth = np.convolve(forces, np.ones((window,))/window, mode='valid')
    return dist_smooth, forces_smooth


def find_derivative(dist, forces):
    forces_derivative = np.diff(np.array(forces))/np.diff(np.array(dist))
    dist_derivative = np.array(dist[:-1]) + np.diff(np.array(dist))/2
    return running_average(dist_derivative, forces_derivative)


def find_ranges(dist, forces, gap=10):
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
                    and len(current_range) == 1:
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


def fit_last_part(dist, forces, coefficients):
    dlast = Variable('dlast')
    flast = Variable('flast')
    llast = Parameter('llast', value=max(dist) + 20, min=max(dist) + 1)
    plast = Parameter('plast', value=1, min=0.1)
    model = Model({flast: plast * (0.25 / ((1 - dlast / llast) ** 2) - 0.25 + dlast / llast)})
    fit = Fit(model, dlast=dist, flast=forces)
    fit_result = fit.execute()
    coefficients['llast'] = (fit_result.value(llast))
    coefficients['plast'] = (fit_result.value(plast))
    return coefficients


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
        return ranges, dist_smooth, forces_smooth, coefficients
    else:
        raise NameError("Unknown linker for analysis. Known models are " + str(available_linkers) + '.')


def fit_curve(ranges, dist, forces, coefficients, temperature=None, elastic_moduli_boudaries=(250, 400)):
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

    coefficients = fit_last_part(dist_to_fit[-1], forces_to_fit[-1], coefficients)

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


def find_rupture_forces(dist, forces, ranges):
    rupture_forces = [max([forces[_] for _ in range(len(dist))
                           if current_range[0] < dist[_] < current_range[1]]) for current_range in ranges]
    return rupture_forces


def find_total_stretching(coefficients, show_plot):
    lengths = [max(coeficient_set['L'][1:]) for coeficient_set in coefficients]
    (mu, sigma) = norm.fit(lengths)
    n, bins, patches = plt.hist(lengths, density=True, facecolor='green', alpha=0.75)
    x = np.linspace(min(lengths), max(lengths), 100)
    plt.plot(x, norm.pdf(x, mu, sigma), 'r--', linewidth=2)

    # plot
    plt.xlabel('Lenghts [nm]')
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Histogram\ of\ protein\ lengths:}\ \mu=%.3f,\ \sigma=%.3f$' % (mu, sigma))
    plt.grid(True)
    if show_plot:
        plt.show()
        plt.close('all')
    else:
        plt.savefig('total_length_histogram.png')
    return round(mu, 3)


def make_contour_gain_histo(coefficients, show_plot):
    to_plot = []
    colors = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta']
    for coefficient_line in coefficients:
        data = coefficient_line['L']
        for n in range(1, len(data)):
            if len(to_plot) < n:
                to_plot.append([])
            if n == 1:
                to_plot[n-1].append(data[n])
            else:
                to_plot[n-1].append(data[n]-data[n-1])
    results = []
    for k in range(len(to_plot)):
        n, bins, patches = plt.hist(to_plot[k], density=True, facecolor=colors[k % len(colors)], alpha=0.5)
        (mu, sigma) = norm.fit(to_plot[k])
        results.append(round(mu, 3))
        x = np.linspace(min(to_plot[k]), max(to_plot[k]), 100)
        plt.plot(x, norm.pdf(x, mu, sigma), colors[k % len(colors)][0] + '--', linewidth=2)
    plt.xlabel('Countour length gain [nm]')
    plt.ylabel('Probability')
    plt.title('Histogram of contour length gain')
    plt.grid(True)
    if show_plot:
        plt.show()
        plt.close('all')
    else:
        plt.savefig('contour_length_gain_histogram.png')
    return results


def make_force_histogram(forces, show_plot, cutoff=0):
    if cutoff > 0:
        forces = [f for f in forces if f <= cutoff]
    n, bins, patches = plt.hist(forces, density=True, facecolor='blue', alpha=0.5)
    plt.xlabel('Rupture forces [pN]')
    plt.ylabel('Probability')
    plt.title('Histogram of rupture forces')
    plt.grid(True)
    if show_plot:
        plt.show()
        plt.close('all')
    else:
        plt.savefig('rupture_forces_histogram.png')
    return n, bins


def transform(d, f, length, pdna, pp, stiff):
    result = (d - length * (1-1/np.sqrt(pdna*f) + f/stiff)) / (1 - 1/np.sqrt(pp*f))
    return result


def transform_coordinates(dist, forces, coefficients):
    result = []
    length = coefficients['L'][0]
    pdna = coefficients['pDNA/KbT'],
    pp = coefficients['pProtein/KbT']
    stiff = coefficients['K']
    for d, f in zip(dist, forces):
        result.append((transform(d, f, length, pdna, pp, stiff),f))
    return result


def plot_transformed_coordinates(contour_lengths, show_plot):
    n, bins, patches = plt.hist(contour_lengths, density=True, facecolor='blue', alpha=0.5)
    plt.xlabel('Contour length [nm]')
    plt.ylabel('Probability')
    plt.title('Histogram of contour lengths')
    plt.grid(True)
    if show_plot:
        plt.show()
        plt.close('all')
    else:
        plt.savefig('contour_length_histogram.png')
    return n, bins


def calculate_protein_force(d, L, p):
    return p * (0.25/((1-d/L)**2) - 0.25 + d/L)


def calculate_dna(F, L, K, p):
    return L * (1 - np.sqrt(p/F) + F/K)


def calculate_protein(F, L, p):
    return L * (1 - np.sqrt(p/F))


def calculate_total(F, Ld, K, pd, Lp, pp):
    return calculate_dna(F, Ld, K, pd) + calculate_protein(F, Lp, pp)


def cluster_coefficients(coefficients, maxgap=15):
    coefficients.sort()
    clusters = [[coefficients[0]]]
    for x in coefficients[1:]:
        if abs(x - clusters[-1][-1]) <= maxgap:
            clusters[-1].append(x)
        else:
            clusters.append([x])
    singles = [k for k in range(len(clusters)) if len(clusters[k]) <= 4]
    clusters.append([loc for k in singles for loc in clusters[k]])
    for k in list(reversed(singles)):
        clusters.pop(k)
    return clusters


def plot_coefficients(dist_total, forces_total, ranges_total, coefficients_total, name, linker='None', show_plots=False,
                      columns=4, cluster_max_gap=15, residues_distance=3.88):
    rows = int(np.ceil(float(len(ranges_total))/columns))
    fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
    colors = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta']
    results = []
    for k in range(len(ranges_total) + 1):
        row = int(k/4)
        col = k % 4
        if k == 0:
            axes[row, col].set_title('Histogram')
            axes[row, col].set_xlabel('Countour length gain [A]')
            axes[row, col].set_ylabel('Probability')
            all_coefficients = [coefficient - r['L'][0] for r in coefficients_total for coefficient in r['L'][1:]]
            coefficients_clusters = cluster_coefficients(all_coefficients, maxgap=cluster_max_gap)
            min_l = min(all_coefficients)
            max_l = max(all_coefficients)
            axes[row, col].set_xlim(min_l, max_l)
            x = np.linspace(min_l, max_l, 100)
            for r in range(len(coefficients_clusters)-1):
                n, bins, patches = axes[row, col].hist(coefficients_clusters[r], density=True,
                                                       facecolor=colors[r % len(colors)], alpha=0.5)
                (mu, sigma) = norm.fit(coefficients_clusters[r])
                residues = int(round(mu/residues_distance, 0))
                results.append([round(mu, 3), residues])
                label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
                axes[row, col].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)][0] + '--', label=label)
            n, bins, patches = axes[row, col].hist(coefficients_clusters[-1], bins=np.arange(min_l, max_l+5, 5),
                                                   density=True, facecolor='k', alpha=0.5)
        else:
            dist = dist_total[k-1]
            forces = forces_total[k-1]
            ranges = ranges_total[k-1]
            coefficients = coefficients_total[k-1]
            axes[row, col].set_xlim(min(dist), max(dist))
            axes[row, col].set_ylim(0, max(forces))
            axes[row, col].set_title(name + '_' + str(k))
            axes[row, col].set_xlabel('Extension [A]')
            axes[row, col].set_ylabel('Force [pN]')
            axes[row, col].plot(np.array(dist), np.array(forces))
            for l in range(len(coefficients['L'])):
                length = coefficients['L'][l]
                residues = int(round(length / residues_distance, 0))
                label = str(round(length, 3)) + ' (' + str(residues) + 'AA)'
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


def save_data(contour_length_gain):
    return