""" The class of whole measurements corresponding to one structure """

from datetime import datetime
from curve_parsing import Curve
from tools import *
import matplotlib.pyplot as plt
from scipy.stats import cauchy, ks_2samp
import numpy as np


class Structure:
    def __init__(self, info, parameters):
        self.info = info
        self.parameters = parameters
        self.curves = []
        self.dists = []
        self.forces = []
        self.coefficients = []
        self.hist_ranges = []
        self.rupture_forces = {}
        self.lengths = []
        self.pprot = 0
        self.contour_length = []
        return

    def _extract_curves(self):
        """ The function extracting the data from file"""
        fname = self.parameters['data_path'] + self.parameters['data_file_prefix'] + \
                self.info['name'] + self.parameters['data_file_suffix']

        with open(fname, 'r') as myfile:
            # reading values
            values = []
            for line in myfile.readlines():
                if 'd' not in line and 'D' not in line and 'TIME' not in line:
                    values.append(line.strip().split(';'))
        # extracting traces
        beg = len(values[0]) % 2
        end = len(values[0])
        for trace in range(beg, end, 2):
            pairs = [[float(values[point][trace]), float(values[point][trace + 1])]
                     for point in range(len(values)) if values[point][trace]]
            pairs = sorted(pairs)
            dist = [v[0] for v in pairs]
            force = [v[1] for v in pairs]
            if self.info['unit'] == 'nm':
                dist = [10 * d for d in dist]
            self.dists.append(dist)
            self.forces.append(force)
        return

    def _make_partial_plots(self, debug=False):
        if debug:
            print("\t-> " + str(datetime.now().time()) + " Making plots of indicidual fits.")
        number = len(self.curves)       # number of traces to plot
        columns = self.parameters['columns']
        rows = max(int(np.ceil(float(number) / columns)), 2)
        fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
        for k in range(number):
            row = int(k/4)
            col = k % 4
            self.curves[k].plot(position=axes[row, col])
        fig.tight_layout()
        if debug:
            print("Saving contour lengths figure to " + self.info['name'] + '_contour_lengths.png')
        plt.savefig(self.info['name'] + '_contour_lengths.png')
        return

    def _gather_contour_lengths(self):
        return sorted([length for c in self.curves for length in c.contourlengths])

    def _gather_rupture_forces(self):
        rupture_forces = {}
        for c in self.curves:
            rupture_forces = merge_dicts(rupture_forces, c.rupture_forces)
        return rupture_forces

    def _plot_contour_length_histo(self, position):
        lengths_range = np.linspace(min(self.contour_length), max(self.contour_length))
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')
        for r in range(len(self.hist_ranges)):
            current_range = self.hist_ranges[r]
            values = [clength for clength in self.contour_length if current_range[0] <= clength <= current_range[1]]
            n, bins, patches = position.hist(values, bins=200, color=colors[r % len(colors)])
            area = find_area(n, bins)
            mu, gamma = self.lengths[r][:2]
            label = 'L=' + str(round(mu, 3))
            position.plot(lengths_range, area * cauchy.pdf(lengths_range, mu, gamma), linestyle='--', label=label,
                          color=colors[r % len(colors)])
        position.legend()
        return

    def _make_histograms(self, debug=False):
        if debug:
            print("\t-> " + str(datetime.now().time()) + " Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        # the contour length histogram
        self._plot_contour_length_histo(axes[0, 0])

        # distance = self.parameters['residues_distance'][self.info['source']]
        # dudko_analysis = {}
        # persistence = sum([c.coefficients['pprot'] for c in self.curves])/len(self.curves)
        #
        # # contour length histogram
        # axes[0, 0].set_title('Contour length histogram')
        # axes[0, 0].set_xlabel('Countour length [A]')
        # axes[0, 0].set_ylabel('Probability')
        # all_coefficients = [coefficient for c in self.curves for coefficient in c.coefficients['L']]
        # coefficients_clusters = cluster_coefficients(all_coefficients, maxgap=self.parameters['cluster_max_gap'])
        # min_l = min(all_coefficients)
        # max_l = max(all_coefficients)
        # axes[0, 0].set_xlim(min_l, max_l)
        # x = np.linspace(min_l, max_l, 100)
        # for r in range(len(coefficients_clusters)-1):
        #     n, bins, patches = axes[0, 0].hist(coefficients_clusters[r], density=True, facecolor=colors[r % len(colors)],
        #                                        alpha=0.5)
        #     (mu, sigma) = norm.fit(coefficients_clusters[r])
        #     residues = int(round(mu / distance, 0)) + 1
        #     self.coefficients.append([round(mu, 3), round(sigma, 3), residues])
        #     label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
        #     axes[0, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
        #     dudko_analysis[r] = {'mu': mu}
        # n, bins, patches = axes[0, 0].hist(coefficients_clusters[-1], bins=np.arange(min_l, max_l + 5, 5),
        #                                        density=True, facecolor='k', alpha=0.5)
        # axes[0, 0].legend()
        # if debug:
        #     print("\t\t-> " + str(datetime.now().time()) + " Done contour length histogram.")
        #
        # # forces histogram
        # axes[1, 0].set_title('Rupture forces histogram')
        # axes[1, 0].set_xlabel('Rupture force [pN]')
        # axes[1, 0].set_ylabel('Probability')
        # rupture_forces = self._gather_rupture_forces()
        # min_f = min(list(rupture_forces.values()))
        # max_f = max(list(rupture_forces.values()))
        # axes[1, 0].set_xlim(min_f, max_f)
        # x = np.linspace(min_f, max_f, 100)
        # for r in range(len(coefficients_clusters)-1):
        #     force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[r]]
        #     n, bins, patches = axes[1, 0].hist(force_cluster, density=True, facecolor=colors[r % len(colors)], alpha=0.5)
        #     dudko_analysis[r]['n'] = n
        #     dudko_analysis[r]['bins'] = bins
        #     (mu, sigma) = norm.fit(force_cluster)
        #     label = str(round(mu, 3)) # + ' (L=' + str(results['contour_length_histo'][r][0]) + ')'
        #     axes[1, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
        # force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[-1]]
        # n, bins, patches = axes[1, 0].hist(force_cluster, bins=np.arange(min_l, max_l + 5, 5),
        #                                        density=True, facecolor='k', alpha=0.5)
        # axes[1, 0].legend()
        # if debug:
        #     print("\t\t-> " + str(datetime.now().time()) + " Done rupture forces histogram.")
        #
        # # transformed contour length histogram
        # axes[0, 1].set_title('Transformed contour length histogram')
        # axes[0, 1].set_xlabel('Contour length [pN]')
        # axes[0, 1].set_ylabel('Probability')
        # contour_lengths = self._gather_contour_lengths()
        # fit_ranges = find_hist_ranges(contour_lengths)
        # min_l = min(contour_lengths)
        # max_l = max(contour_lengths)
        # axes[0, 1].set_xlim(min_l, max_l)
        # x = np.linspace(min_l, max_l, 100)
        # for r in range(len(fit_ranges)):
        #     data_to_fit = [value for value in contour_lengths if fit_ranges[r][0] <= value <= fit_ranges[r][1]]
        #     n, bins, patches = axes[0, 1].hist(data_to_fit, density=True, facecolor=colors[r % len(colors)], alpha=0.5)
        #     (mu, sigma) = norm.fit(data_to_fit)
        #     residues = int(round(mu / distance, 0)) + 1
        #     self.lengths.append([round(mu, 3), round(sigma, 3), residues])
        #     label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
        #     axes[0, 1].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
        # axes[0, 1].legend()
        # if debug:
        #     print("\t\t-> " + str(datetime.now().time()) + " Done histogram of transformed coordinates.")
        #
        # # Dudko analysis
        # axes[1, 1].set_title('Dudko analysis')
        # axes[1, 1].set_xlabel('Rupture force [pN]')
        # axes[1, 1].set_ylabel('Rupture time [s]')
        # axes[1, 1].set_xlim(min_f, max_f)
        # dudko_result = []
        # for r in dudko_analysis.keys():
        #     l = dudko_analysis[r]['mu']
        #     heights = np.array(dudko_analysis[r]['n'])
        #     heights_div = np.array([heights[x] for x in range(len(heights)) if heights[x] != 0])
        #     bins = dudko_analysis[r]['bins']
        #     df = bins[1] - bins[0]
        #     plt_range = [k for k in range(1, len(bins)) if heights[k - 1] != 0]
        #     x = np.array([bins[0] + (k - 1 / 2) * df for k in plt_range])
        #     load = force_load(x, persistence, l, self.parameters['speed'])
        #     area = np.array(
        #         [heights[k - 1] / 2 + sum([heights[i - 1] for i in range(k + 1, len(bins))]) for k in plt_range])
        #     y = (area * df) / (load * heights_div)
        #     label = 'L=' + str(round(l, 3))
        #     dudko_result.append(fit_dudko(x, y))
        #     x_smooth = np.linspace(min(x), max(x))
        #     t0 = dudko_result[-1]['1/2']['t0']
        #     gamma = dudko_result[-1]['1/2']['x'] / (2 * dudko_result[-1]['1/2']['g'])
        #     g = dudko_result[-1]['1/2']['g']
        #     y_smooth = t0 * (1 - gamma * x_smooth) ** (-1) * np.exp(-g * (1 - (1 - gamma * x_smooth) ** 2))
        #     axes[1, 1].plot(x, y, 'o', color=colors[r % len(colors)], label=label)
        #     axes[1, 1].plot(x_smooth, y_smooth, '-', color=colors[r % len(colors)])
        # axes[1, 1].legend()
        fig.tight_layout()
        if debug:
            # print("\t\t-> " + str(datetime.now().time()) + " Done analysis of Dudko equation.")
            print("Saving histograms figure to " + self.info['name'] + '_histograms.png')
        plt.savefig(self.info['name'] + '_histograms.png')
        return

    def make_plots(self, debug=False):
        self._make_partial_plots(debug=debug)
        self._make_histograms(debug=debug)
        return

    def save_data(self):
        result = []
        separator = '################\n'

        # general info
        result.append('General info:')
        result.append('Name:\t\t\t' + self.info['name'])
        result.append('Residues:\t\t\t' + str(self.info['residues']))
        result.append('End-to-end distance:\t' + str(self.info['distance']))
        result.append('Linker:\t\t\t' + str(self.info['linker']))
        result.append('Unit:\t\t\t' + str(self.info['unit']))
        result.append('Data source:\t\t' + str(self.info['source']))
        result.append('Pulling speed:\t\t\t' + str(self.info['speed']))
        result.append(separator)

        # parameters
        result.append('Calculation parameters:')
        result.append('Data path:\t\t' + str(self.parameters['data_path']))
        result.append('Data file prefix:\t\t' + str(self.parameters['data_file_prefix']))
        result.append('Data file suffix:\t\t' + str(self.parameters['data_file_suffix']))
        result.append('Residue-residue distance:\t' + str(self.parameters['residues_distance']))
        result.append('Minimal distance between jumps:\t\t' + str(self.parameters['minimal_stretch_distance']))
        result.append('Low force cutoff:\t\t' + str(self.parameters['low_force_cutoff']))
        result.append('High force cutoff:\t\t' + str(self.parameters['high_force_cutoff']))
        result.append('Minimal gap between peaks in cluster:\t\t' + str(self.parameters['cluster_max_gap']))
        result.append('Number of columns in individual plots:\t\t' + str(self.parameters['columns']))
        result.append(separator)

        # summary of individual curve
        result.append('Summary of individual curves')
        for k in range(len(self.curves)):
            result.append(str(k) + '/' + str(len(self.curves)))
            result.append(self.curves[k].summary())
        result.append(separator)

        # summary of the cumulative statistics
        result.append('Summary of the cummulative statistics')
        result.append('pProt:\t\t' + str(self.pprot))
        result.append('Contour length\tgamma\tks pValue')
        for mu, gamma, pvalue in self.lengths:
            result.append(str(mu) + '\t\t' + str(gamma) + '\t' + str(pvalue))
        result.append('Contour length histogram delimiting regions:')
        result.append(str(self.hist_ranges))



        fname = self.info['name'] + '_results'
        with open(fname, 'w') as myfile:
            myfile.write('\n'.join(result))
        return

    def _calcualate_contour_lengths(self, pprot):
        self.contour_length = []
        for k in range(len(self.dists)):
            dist = np.array([self.dists[k][_] for _ in range(len(self.dists[k]))
                             if self.forces[k][_] > self.parameters['low_force_cutoff']])
            xs = np.array([invert_wlc_np(f, pprot) for f in self.forces[k] if
                           f > self.parameters['low_force_cutoff']])
            self.contour_length += list(dist / xs)
        self.contour_length = np.array(self.contour_length)
        return

    def _find_constants(self):
        pprots = [c.coefficients['pprot'] for c in self.curves]
        pprot = sum(pprots)/len(pprots)
        self._calcualate_contour_lengths(pprot)
        self.hist_ranges = find_hist_ranges(self.contour_length)
        for r in range(len(self.hist_ranges)):
            current_range = self.hist_ranges[r]
            values = [clength for clength in self.contour_length if current_range[0] <= clength <= current_range[1]]
            mu, gamma = cauchy.fit(values)
            test_statistic = cauchy.rvs(mu, gamma, len(values))
            result = [mu, gamma, ks_2samp(values, test_statistic).pvalue]
            self.lengths.append(result)
        self.pprot = pprot
        return

    def analyze(self, debug=False):
        """ The main function for the analysis of particular structure."""
        self._extract_curves()
        if debug:
            print(str(datetime.now().time()) + "\tAnalyzing " + self.info['name'] +
                  ". Found " + str(len(self.dists)) + " traces to analyze. Analyzing:")
        for k in range(5):  #range(len(self.dists)):
            if debug:
                print(str(datetime.now().time()) + "\t" + str(len(self.curves) + 1) + "/" + str(len(self.dists)))
            self.curves.append(Curve(self.dists[k], self.forces[k], self.info, self.parameters, debug=debug))
            self.curves[-1].analyze()
        self._find_constants()
        self.make_plots(debug=debug)
        self.save_data()
        return
