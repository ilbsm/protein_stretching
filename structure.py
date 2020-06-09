""" The class of whole measurements corresponding to one structure """

from datetime import datetime
from curve_parsing import Curve
from tools import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import norm
import numpy as np

colors = [mcolors.CSS4_COLORS['red'],
          mcolors.CSS4_COLORS['green'],
          mcolors.CSS4_COLORS['blue'],
          mcolors.CSS4_COLORS['yellow'],
          mcolors.CSS4_COLORS['cyan'],
          mcolors.CSS4_COLORS['magenta'],
          mcolors.CSS4_COLORS['orange'],
          mcolors.CSS4_COLORS['purple'],
          mcolors.CSS4_COLORS['lime']]


class Structure:
    def __init__(self, info, parameters):
        self.info = info
        self.parameters = parameters
        self.curves = []
        self.dists = []
        self.forces = []
        self.coefficients = []
        self.rupture_forces = {}
        self.lengths = []
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
        number = len(self.curves)       # number of traces to plot
        columns = self.parameters['columns']
        distance = self.parameters['residues_distance'][self.info['source']]
        rows = max(int(np.ceil(float(number) / columns)), 2)
        fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
        for k in range(number):
            row = int(k/4)
            col = k % 4
            dist = self.dists[k]
            forces = self.forces[k]
            ranges = self.curves[k].ranges
            coefficients = self.curves[k].coefficients
            pprot = coefficients['pprot']
            f_plot = np.linspace(max(min(forces), 0.01), max(forces))

            axes[row, col].set_xlim(min(dist), max(dist))
            axes[row, col].set_ylim(0, max(forces))
            axes[row, col].set_title(self.info['name'] + '_' + str(k+1))
            axes[row, col].set_xlabel('Extension [A]')
            axes[row, col].set_ylabel('Force [pN]')
            axes[row, col].plot(np.array(dist), np.array(forces))

            for l in range(len(coefficients['L'])):
                length = coefficients['L'][l]
                residues = int(round(length / distance, 0)) + 1
                label = str(round(length, 3)) + ' (' + str(residues) + ' AA)'
                d_plot = length * np.array([invert_wlc_np(f, pprot) for f in f_plot])
                if 'ldna' in coefficients.keys():
                    ldna = coefficients['ldna']
                    pdna = coefficients['pdna']
                    ddna = ldna * np.array([invert_wlc_np(f, pdna) for f in f_plot])
                    d_plot += 2*ddna
                axes[row, col].plot(d_plot, f_plot, label=label, color=colors[l % len(colors)])
            for r in ranges:
                axes[row, col].axvline(x=r[0], linestyle='--', color='#808080', linewidth=0.5)
                axes[row, col].axvline(x=r[1], linestyle='--', color='#808080', linewidth=0.5)
            axes[row, col].legend()
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

    def _make_histograms(self, debug=False):
        if debug:
            print("\t-> " + str(datetime.now().time()) + " Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))
        distance = self.parameters['residues_distance'][self.info['source']]
        dudko_analysis = {}
        persistence = sum([c.coefficients['pprot'] for c in self.curves])/len(self.curves)

        # contour length histogram
        axes[0, 0].set_title('Contour length histogram')
        axes[0, 0].set_xlabel('Countour length [A]')
        axes[0, 0].set_ylabel('Probability')
        all_coefficients = [coefficient for c in self.curves for coefficient in c.coefficients['L']]
        coefficients_clusters = cluster_coefficients(all_coefficients, maxgap=self.parameters['cluster_max_gap'])
        min_l = min(all_coefficients)
        max_l = max(all_coefficients)
        axes[0, 0].set_xlim(min_l, max_l)
        x = np.linspace(min_l, max_l, 100)
        for r in range(len(coefficients_clusters)-1):
            n, bins, patches = axes[0, 0].hist(coefficients_clusters[r], density=True, facecolor=colors[r % len(colors)],
                                               alpha=0.5)
            (mu, sigma) = norm.fit(coefficients_clusters[r])
            residues = int(round(mu / distance, 0)) + 1
            self.coefficients.append([round(mu, 3), round(sigma, 3), residues])
            label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
            axes[0, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
            dudko_analysis[r] = {'mu': mu}
        n, bins, patches = axes[0, 0].hist(coefficients_clusters[-1], bins=np.arange(min_l, max_l + 5, 5),
                                               density=True, facecolor='k', alpha=0.5)
        axes[0, 0].legend()
        if debug:
            print("\t\t-> " + str(datetime.now().time()) + " Done contour length histogram.")

        # forces histogram
        axes[1, 0].set_title('Rupture forces histogram')
        axes[1, 0].set_xlabel('Rupture force [pN]')
        axes[1, 0].set_ylabel('Probability')
        rupture_forces = self._gather_rupture_forces()
        min_f = min(list(rupture_forces.values()))
        max_f = max(list(rupture_forces.values()))
        axes[1, 0].set_xlim(min_f, max_f)
        x = np.linspace(min_f, max_f, 100)
        for r in range(len(coefficients_clusters)-1):
            force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[r]]
            n, bins, patches = axes[1, 0].hist(force_cluster, density=True, facecolor=colors[r % len(colors)], alpha=0.5)
            dudko_analysis[r]['n'] = n
            dudko_analysis[r]['bins'] = bins
            (mu, sigma) = norm.fit(force_cluster)
            label = str(round(mu, 3)) # + ' (L=' + str(results['contour_length_histo'][r][0]) + ')'
            axes[1, 0].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
        force_cluster = [rupture_forces[coefficient] for coefficient in coefficients_clusters[-1]]
        n, bins, patches = axes[1, 0].hist(force_cluster, bins=np.arange(min_l, max_l + 5, 5),
                                               density=True, facecolor='k', alpha=0.5)
        axes[1, 0].legend()
        if debug:
            print("\t\t-> " + str(datetime.now().time()) + " Done rupture forces histogram.")

        # transformed contour length histogram
        axes[0, 1].set_title('Transformed contour length histogram')
        axes[0, 1].set_xlabel('Contour length [pN]')
        axes[0, 1].set_ylabel('Probability')
        contour_lengths = self._gather_contour_lengths()
        fit_ranges = find_hist_ranges(contour_lengths)
        min_l = min(contour_lengths)
        max_l = max(contour_lengths)
        axes[0, 1].set_xlim(min_l, max_l)
        x = np.linspace(min_l, max_l, 100)
        for r in range(len(fit_ranges)):
            data_to_fit = [value for value in contour_lengths if fit_ranges[r][0] <= value <= fit_ranges[r][1]]
            n, bins, patches = axes[0, 1].hist(data_to_fit, density=True, facecolor=colors[r % len(colors)], alpha=0.5)
            (mu, sigma) = norm.fit(data_to_fit)
            residues = int(round(mu / distance, 0)) + 1
            self.lengths.append([round(mu, 3), round(sigma, 3), residues])
            label = str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
            axes[0, 1].plot(x, norm.pdf(x, mu, sigma), colors[r % len(colors)], linestyle='--', label=label)
        axes[0, 1].legend()
        if debug:
            print("\t\t-> " + str(datetime.now().time()) + " Done histogram of transformed coordinates.")

        # Dudko analysis
        axes[1, 1].set_title('Dudko analysis')
        axes[1, 1].set_xlabel('Rupture force [pN]')
        axes[1, 1].set_ylabel('Rupture time [s]')
        axes[1, 1].set_xlim(min_f, max_f)
        dudko_result = []
        for r in dudko_analysis.keys():
            l = dudko_analysis[r]['mu']
            heights = np.array(dudko_analysis[r]['n'])
            heights_div = np.array([heights[x] for x in range(len(heights)) if heights[x] != 0])
            bins = dudko_analysis[r]['bins']
            df = bins[1] - bins[0]
            plt_range = [k for k in range(1, len(bins)) if heights[k - 1] != 0]
            x = np.array([bins[0] + (k - 1 / 2) * df for k in plt_range])
            load = force_load(x, persistence, l, self.parameters['speed'])
            area = np.array(
                [heights[k - 1] / 2 + sum([heights[i - 1] for i in range(k + 1, len(bins))]) for k in plt_range])
            y = (area * df) / (load * heights_div)
            label = 'L=' + str(round(l, 3))
            dudko_result.append(fit_dudko(x, y))
            x_smooth = np.linspace(min(x), max(x))
            t0 = dudko_result[-1]['1/2']['t0']
            gamma = dudko_result[-1]['1/2']['x'] / (2 * dudko_result[-1]['1/2']['g'])
            g = dudko_result[-1]['1/2']['g']
            y_smooth = t0 * (1 - gamma * x_smooth) ** (-1) * np.exp(-g * (1 - (1 - gamma * x_smooth) ** 2))
            axes[1, 1].plot(x, y, 'o', color=colors[r % len(colors)], label=label)
            axes[1, 1].plot(x_smooth, y_smooth, '-', color=colors[r % len(colors)])
        axes[1, 1].legend()
        fig.tight_layout()
        if debug:
            print("\t\t-> " + str(datetime.now().time()) + " Done analysis of Dudko equation.")
            print("Saving histograms figure to " + self.info['name'] + '_histograms.png')
        plt.savefig(self.info['name'] + '_histograms.png')
        return

    def make_plots(self, debug=False):
        self._make_partial_plots(debug=debug)
        self._make_histograms(debug=debug)
        return

    def save_data(self):
        return

    def analyze(self, debug=False):
        """ The main function for the analysis of particular structure."""
        self._extract_curves()
        if debug:
            print(str(datetime.now().time()) + "\tAnalyzing " + self.info['name'] +
                  ". Found " + str(len(self.dists)) + " traces to analyze. Analyzing:")
        for k in range(len(self.dists)):
            if debug:
                print(str(datetime.now().time()) + "\t" + str(len(self.curves) + 1) + "/" + str(len(self.dists)))
            self.curves.append(Curve(self.dists[k], self.forces[k], self.info, self.parameters, debug=debug))
            self.curves[-1].fit()
            self.curves[-1].transform_coordinates()
            self.curves[-1].find_rupture_forces()
            # self.curves[-1].find_energies()
        self.make_plots(debug=debug)
        # self.save_data()
        return


