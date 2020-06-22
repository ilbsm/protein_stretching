from datetime import datetime
from scipy.stats import cauchy, ks_2samp
from .stretching_tools import *
import pandas as pd


class Trace:
    def __init__(self, number, data, parameters=None, debug=False):
        headers = {'d_' + str(number): 'd', 'F_' + str(number): 'F'}
        self.name = number
        self.parameters = parameters
        self.debug = debug

        self.data = data.rename(columns=headers)
        self.data = self.data[self.data['F'] > self.parameters['low_force_cutoff']].sort_values(by='d')

        self.smoothed = pd.DataFrame()
        self.derivatives = pd.DataFrame()

        self.ranges = []
        self.hist_ranges = []
        self.max_contour_length = 0
        self.coefficients = {}
        self.rupture_forces = {}
        return

    def _generate_smooths(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Generating smooth trace...")
        f_smooth = []
        d_smooth = []

        # averaging over the interval
        d_range = np.linspace(min(self.data['d']), max(self.data['d']), 101)
        for k in range(len(d_range)-1):
            d_min = d_range[k]
            d_max = d_range[k+1]
            forces_interval = self.data[self.data['d'].between(d_min, d_max, inclusive=True)]['F']
            if len(forces_interval) > 0:
                d_smooth.append(d_min)
                f_smooth.append(sum(forces_interval) / len(forces_interval))

        d_smooth, f_smooth = running_average(d_smooth, f_smooth)
        self.smoothed = pd.DataFrame({'d': d_smooth, 'F': f_smooth})
        return

    def _generate_derivative(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Generating trace derivatives...")
        d_diff, f_diff = find_derivative(self.smoothed['d'], self.smoothed['F'])
        self.derivatives = pd.DataFrame({'d': d_diff, 'F': f_diff})
        return

    def _generate_ranges(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Generating trace ranges...")
        ranges = []
        current_range = [self.derivatives['d'][0]]
        neg_diff = -1
        force_cutoff = self.parameters['high_force_cutoff']
        for d in range(len(self.derivatives['d'])):
            if self.derivatives['F'][d] < 0:
                if neg_diff < 0:
                    neg_diff = d
                # if the window of negative derivative is at least equal 'minimal_stretch_distance' parameter,
                # we finish the current range
                if neg_diff > 0 and self.derivatives['d'][d] - self.derivatives['d'][neg_diff] >= \
                        self.parameters['minimal_stretch_distance'] and len(current_range) == 1 \
                        and self.smoothed['F'][d] >= force_cutoff:
                    current_range.append(self.derivatives['d'][neg_diff])
                    ranges.append(current_range)
                    current_range = []
            # establishing new range beginning as the first point when derivative becomes again positive
            if self.derivatives['F'][d] > 0:
                if neg_diff > 0 and len(current_range) == 0:
                    current_range.append(self.derivatives['d'][d])
                neg_diff = -1
        if len(current_range) == 1:
            current_range.append(self.derivatives.tail(1)['d'].values[0])
            ranges.append(current_range)
        self.ranges = ranges
        return

    def _fit_cl_dna_experiment(self):
        pprot = self.parameters['initial_guess']['pprot']
        pdna = self.parameters['initial_guess']['pdna']
        elast = self.parameters['initial_guess']['K']
        ldna = self.parameters['initial_guess']['ldna']

        self.data['x_dna'] = np.array([invert_wlc(f, pdna, elast) for f in self.data['F']])
        self.data['x_prot'] = np.array([invert_wlc(f, pprot) for f in self.data['F']])
        self.data['d_dna'] = ldna * self.data['x_dna']
        self.data['d_prot'] = self.data['d'] - self.data['d_dna']
        self.data['hist_values'] = self.data['d_prot'] / self.data['x_prot']
        self.hist_ranges = find_hist_ranges(self.data['hist_values'])
        self.max_contour_length = max(self.data['hist_values'])

        parameters = []
        for current_range in self.hist_ranges:
            values = [v for v in self.data['hist_values'] if current_range[0] <= v <= current_range[1]]
            mu, gamma = cauchy.fit(values)
            test_statistic = cauchy.rvs(mu, gamma, len(values))
            parameters.append([mu, gamma, ks_2samp(values, test_statistic).pvalue])
        coefficients = {'pprot': pprot, 'pdna': pdna, 'ldna': ldna, 'K': elast, 'lengths': parameters}
        return coefficients

    def _fit_cl_dna_theory(self):
        coefficients = {}
        return coefficients

    def _fit_cl_none_experiment(self):
        coefficients = {}
        return coefficients

    def _fit_cl_none_theory(self):
        coefficients = {}
        return coefficients

    def fit_contour_lengths(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Fitting contour lengths...")
        if self.parameters['linker'] == 'dna' and self.parameters['source'] == 'experiment':
            coefficients = self._fit_cl_dna_experiment()
        elif self.parameters['linker'] == 'dna' and self.parameters['source'] == 'theory':
            coefficients = self._fit_cl_dna_theory()
        elif self.parameters['linker'] == 'none' and self.parameters['source'] == 'experiment':
            coefficients = self._fit_cl_none_experiment()
        elif self.parameters['linker'] == 'none' and self.parameters['source'] == 'theory':
            coefficients = self._fit_cl_none_theory()
        else:
            raise ValueError("Unknown combination of data source and linker. \n"
                             "Got data source " + self.parameters['source'] + " and linker " +
                             self.parameters['linker'])
        self.coefficients = coefficients
        return

    def find_rupture_forces(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Finding rupture forces...")
        # for r in range(len(self.ranges)):
        #     current_range = self.ranges[r]
        #     length = self.coefficients['L'][r]
        #     try:
        #         force = max([self.data['F'][_] for _ in range(len(self.data['d']))
        #                      if current_range[0] <= self.data['d'][_] <= current_range[1]])
        #     except ValueError:
        #         force = -1
        #     self.rupture_forces[force] = length
        return

    def find_energies(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + "\t Finding energies...")
        return

    # Plotting #
    def plot_histogram(self, position=None, max_contour_length=None):
        position.set_title('Contour length histogram (trace ' + str(self.name) + ')')
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Counts')
        if max_contour_length:
            position.set_xlim(0, max_contour_length)
        lspace = np.linspace(self.hist_ranges[0][0], self.hist_ranges[-1][1])

        for r in range(len(self.hist_ranges)):
            current_range = self.hist_ranges[r]
            mu = self.coefficients['lengths'][r][0]
            gamma = self.coefficients['lengths'][r][1]
            values = [v for v in self.data['hist_values'] if current_range[0] <= v <= current_range[1]]
            residues = int(mu / 0.365)
            label = "L= " + str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
            n, bins, patches = position.hist(values, bins=200)
            area = find_area(n, bins)
            position.plot(lspace, area * cauchy.pdf(lspace, mu, gamma), linestyle='--', label=label)
        position.legend()
        return

    def plot_fits(self, position=None):
        position.set_title('Trace fits (trace ' + str(self.name) + ')')
        position.set_xlim(min(self.data['d']), max(self.data['d']))
        position.set_ylim(0, max(self.data['F']))
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Force [pN]')
        position.plot(self.data['d'], self.data['F'], label='Original data')
        position.plot(self.smoothed['d'], self.smoothed['F'], label='Smoothed data')

        pprot = self.coefficients['pprot']
        pdna = self.coefficients['pdna']
        elast = self.coefficients['K']
        ldna = self.coefficients['ldna']

        for r in range(len(self.ranges)):
            x0, x1 = self.ranges[r]
            position.hlines(max(self.data['F']) / 2, x0 + 2, x1 - 2, color=colors[r % len(colors)], alpha=0.5)
            position.axvline(x=x0, linestyle='--', color='#808080', linewidth=0.5)
            position.axvline(x=x1, linestyle='--', color='#808080', linewidth=0.5)

        f_plot = np.linspace(0.1, max(self.data['F']))
        d_dna = ldna * np.array([invert_wlc(f, pdna, elast) for f in f_plot])
        xprot = np.array([invert_wlc(f, pprot) for f in f_plot])

        for len_nr in range(len(self.coefficients['lengths'])):
            lprot = self.coefficients['lengths'][len_nr][0]
            residues = int(lprot / 0.365)
            d_prot = lprot * xprot
            d_plot = d_dna + d_prot
            label = 'Fit L=' + str(round(lprot, 3)) + ' (' + str(residues) + ' AA)'

            position.plot(d_plot, f_plot, label=label)

        position.legend()
        return

    def analyze(self):
        if self.debug:
            print(str(datetime.now().time()) + "\t Analyzing trace " + str(self.name))
        self._generate_smooths()
        self._generate_derivative()
        self._generate_ranges()
        self.fit_contour_lengths()
        self.find_rupture_forces()
        self.find_energies()
        if self.debug:
            print(str(datetime.now().time()) + "\t Finished trace " + str(self.name))
        return

    def summary(self):
        return
