""" The class of whole measurements corresponding to one structure """

from datetime import datetime
from curve_parsing import Curve
from tools import *
import matplotlib.pyplot as plt
from scipy.stats import cauchy, ks_2samp, norm
import numpy as np


class Structure:
    def __init__(self, info, parameters, cases=None):
        self.info = info
        self.parameters = parameters
        self.curves = []
        self.dists = []
        self.forces = []
        self.coefficients = []
        self.hist_ranges = []
        self.rupture_forces = {}
        self.lengths = {}
        self.pprot = 0
        self.contour_length = []
        self.minf = 0
        self.maxf = 0
        self.dudko_parameters = {}

        self._extract_curves(cases)
        if self.parameters['debug']:
            print(str(datetime.now().time()) + "\tAnalyzing " + self.info['name'] +
                  ". Found " + str(len(self.dists)) + " traces to analyze. Analyzing:")

        for k in range(len(self.dists)):
            if self.parameters['debug']:
                print(str(datetime.now().time()) + "\t" + str(len(self.curves) + 1) + "/" + str(len(self.dists)))
            self.curves.append(Curve(k, self.dists[k], self.forces[k], self.info, self.parameters,
                                     debug=self.parameters['debug']))
        return

    def _extract_curves(self, cases):
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
            if cases and (trace - beg)/2 not in cases:
                continue
            pairs = [[float(values[point][trace]), float(values[point][trace + 1])]
                     for point in range(len(values)) if values[point][trace]]
            pairs = sorted(pairs)
            dist = [v[0] for v in pairs if v[1] > self.parameters['low_force_cutoff'][self.info['source']]]
            force = [v[1] for v in pairs if v[1] > self.parameters['low_force_cutoff'][self.info['source']]]
            # if self.info['unit'] == 'nm':
            #     dist = [10 * d for d in dist]
            self.dists.append(dist)
            self.forces.append(force)
        return

    def _make_partial_plots(self, debug=False):
        if debug:
            print("\t-> " + str(datetime.now().time()) + " Making plots of indicidual fits.")
        number = len(self.curves)       # number of traces to plot
        columns = self.parameters['columns']
        rows = max(int(np.ceil(float(2*number) / columns)), 2)
        fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
        for k in range(0, number, 2):
            self.curves[k].plot_histogram(position=axes[int(k/4), k % 4])
            self.curves[k].plot_fits(position=axes[int((k+1)/4), (k+1) % 4])
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

    def _plot_individual_contour_length_histo(self, position):
        return

    def _plot_total_contour_length_histo(self, position):
        lengths_range = np.linspace(min(self.contour_length), max(self.contour_length))
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')
        for mu in self.lengths.keys():
            current_range = self.lengths[mu]['ranges']
            gamma = self.lengths[mu]['gamma']
            values = [clength for clength in self.contour_length if current_range[0] <= clength <= current_range[1]]
            n, bins, patches = position.hist(values, bins=200, color=self.lengths[mu]['color'])
            area = find_area(n, bins)
            label = 'L=' + str(round(mu, 3))
            position.plot(lengths_range, area * cauchy.pdf(lengths_range, mu, gamma), linestyle='--', label=label,
                          color=self.lengths[mu]['color'])
        position.legend()
        return

    def _plot_forces_histogram(self, position):
        position.set_xlabel('Force [pN]')
        position.set_ylabel('Occurences')
        position.set_title('Rupture force histogram')
        position.set_xlim(self.minf, self.maxf)
        fspace = np.linspace(self.minf, self.maxf, 100)
        for mu in self.lengths.keys():
            forces = self.lengths[mu]['forces']
            n, bins, patches = position.hist(forces, bins=len(forces)/4, density=True,
                                             facecolor=self.lengths[mu]['color'], alpha=0.5)
            mu, sigma = norm.fit(forces)
            label = str(round(mu, 3))
            position.plot(fspace, norm.pdf(fspace, mu, sigma), self.lengths[mu]['color'], linestyle='--', label=label)
            self.lengths[mu]['hist_vals'] = n
            self.lengths[mu]['hist_width'] = bins[1] - bins[0]
            self.lengths[mu]['hist_beg'] = bins[0]
        return

    def _plot_dudko_analysis(self, position):
        position.set_title('Dudko analysis')
        position.set_xlabel('Rupture force [pN]')
        position.set_ylabel('Rupture time [s]')
        position.set_xlim(self.minf, self.maxf)
        for mu in self.lengths:
            self.dudko_parameters[mu] = []
            heights = np.array(self.lengths[mu]['hist_vals'])
            df = self.lengths[mu]['hist_width']
            hist_beg = self.lengths[mu]['hist_beg']
            heights_div = np.array([heights[k] for k in range(len(heights)) if heights[k] != 0])
            x = np.array([hist_beg + (k - 1 / 2) * df for k in range(len(heights)) if heights[k] != 0])
            x_smooth = np.linspace(min(x), max(x))
            load = force_load(x, self.pprot, mu, self.parameters['speed'])
            area = np.array([heights[k - 1] / 2 + sum([heights[i - 1] for i in range(k + 1, len(heights)+1)])
                             for k in range(len(heights)) if heights[k] != 0])
            y = (area * df) / (load * heights_div)
            label = 'L=' + str(round(mu, 3))
            for v in self.parameters['Dudko parameters']:
                self.dudko_parameters[mu][v] = fit_dudko(x, y, v)
                t0 = self.dudko_parameters[mu][v]['t0']
                gamma = self.dudko_parameters[mu][v]['x'] / (2 * self.dudko_parameters[mu][v]['g'])
                g = self.dudko_parameters[mu][v]['g']
                y_smooth = t0 * (1 - gamma * x_smooth) ** (-1) * np.exp(-g * (1 - (1 - gamma * x_smooth) ** 2))
                position.plot(x, y, 'o', color=self.lengths[mu]['color'], label=label)
                position.plot(x_smooth, y_smooth, '-', color=self.lengths[mu]['color'])
        position.legend()
        return

    def _plot_dudko_dependence(self, position):
        position.set_title('Dependence of dG on parameters')
        position.set_xlabel('v')
        position.set_ylabel('dG')
        x = np.array(self.parameters['Dudko parameters'])
        for mu in self.lengths.keys():
            y = np.array([self.dudko_parameters[mu][_]['g'] for _ in x])
            label = str(round(mu, 3))
            position.plot(x, y, color=self.lengths['mu']['color'], label=label)
        position.legend()
        return

    def _make_histograms(self, debug=False):
        if debug:
            print("\t-> " + str(datetime.now().time()) + " Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        # the histogram of contour plots from individual traces
        self._plot_individual_contour_length_histo(axes[0, 0])

        # the total contour length histogram
        self._plot_total_contour_length_histo(axes[0, 1])

        # the rupture forces histogram
        self._plot_forces_histogram(axes[1, 0])

        # Dudko analysis
        self._plot_dudko_analysis(axes[1, 1])

        # # Dudko parameter dependence
        # self._plot_dudko_dependence(axes[1, 1])

        fig.tight_layout()
        if debug:
            # print("\t\t-> " + str(datetime.now().time()) + " Done analysis of Dudko equation.")
            print("Saving histograms figure to " + self.info['name'] + '_histograms.png')
        plt.savefig(self.info['name'] + '_histograms.png')
        return

    def make_plots(self):
        self._make_partial_plots(debug=True)
        # self._make_histograms()
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

    def _find_close_forces(self, r):
        forces = {}
        result = []
        current_range = self.hist_ranges[r]
        for c in self.curves:
            forces = merge_dicts(forces, c.rupture_forces)
        for force in forces.keys():
            if current_range[0] <= forces[force] <= current_range[1]:
                result.append(force)
        return result

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
            self.lengths[mu] = {'gamma': gamma, 'pvalue':  ks_2samp(values, test_statistic).pvalue,
                                'forces': self._find_close_forces(r), 'ranges': current_range,
                                'color': colors[r % len(colors)]}
        self.minf = min([force for mu in self.lengths.keys() for force in self.lengths[mu]['forces']])
        self.maxf = max([force for mu in self.lengths.keys() for force in self.lengths[mu]['forces']])
        self.pprot = pprot
        return

    def analyze(self):
        """ The main function for the analysis of particular structure."""
        for c in self.curves:
            c.analyze()
        # self._find_constants()
        self.make_plots()
        # self.save_data()
        return
