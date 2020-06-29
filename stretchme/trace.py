from datetime import datetime
from scipy.stats import cauchy, ks_2samp, norm
from .stretching_tools import *
import pandas as pd
from sklearn.neighbors import KernelDensity
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema


class Trace:
    def __init__(self, name, data, orig_input, logger=None, parameters=None):
        self.name = name
        self.orig_input = orig_input
        self.logger = logger

        # setting the trace parameters
        self.parameters = parameters

        # filling the data
        self.data = data
        self.data = self.data[self.data['F'] > self.parameters['low_force_cutoff']].sort_values(by='d')
        self.logger.info("Trace name " + str(name) + " got and parsed the data.")

        self.smoothed = pd.DataFrame()

        self.ranges = []
        self.coefficients = {}
        self.rupture_forces = {}
        return

    # generating
    def _generate_smooths(self):
        if self.logger:
            self.logger.info("Generating smooth trace...")
        f_smooth = []
        d_smooth = []

        # averaging over the interval
        d_range = np.linspace(min(self.data['d']), max(self.data['d']), 1001)
        for k in range(len(d_range)-1):
            d_min = d_range[k]
            d_max = d_range[k+1]
            forces_interval = self.data[self.data['d'].between(d_min, d_max, inclusive=True)]['F']
            if len(forces_interval) > 0:
                d_smooth.append(d_min)
                f_smooth.append(forces_interval.mean())

        # running average
        d_smooth, f_smooth = running_average(d_smooth, f_smooth)
        self.smoothed = pd.DataFrame({'d': d_smooth, 'F': f_smooth})
        return

    def _parse_cl_histogram(self, boundary_value=0.001):
        # finding the number of components
        x = np.expand_dims(self.data['hist_values'], 1)
        kde = KernelDensity().fit(x)
        estimator = np.linspace(min(self.data['hist_values']), max(self.data['hist_values']), 1001)
        kde_est = np.exp(kde.score_samples(estimator.reshape(-1, 1)))
        maximas = [estimator[_] for _ in argrelextrema(kde_est, np.greater)[0] if kde_est[_] > boundary_value]

        # finding the range
        # TODO clean it up
        min_boundary = max([estimator[_] for _ in range(len(estimator)) if kde_est[_] < boundary_value and estimator[_] < maximas[0]])
        max_boundary = min([estimator[_] for _ in range(len(estimator)) if kde_est[_] < boundary_value and estimator[_] > maximas[-1]])

        # fitting Gaussians
        x = np.expand_dims(self.data['hist_values'][self.data['hist_values'].between(min_boundary, max_boundary)], 1)
        gmm = GaussianMixture(n_components=len(maximas))  # gmm for two components
        gmm.fit(x)

        # finding parameters
        means = [x[0] for x in gmm.means_]
        widths = [x[0][0] for x in gmm.covariances_]
        heights = [np.exp(gmm.score_samples(np.array(u).reshape(-1, 1))) for u in means]
        boundaries = min_boundary + maximas + max_boundary

        self.ranges = [[boundaries[k], boundaries[k + 1]] for k in range(len(boundaries) - 1)]
        self.coefficients['l_prot'] = [{'mean': means[k], 'width': widths[k], 'height': heights[k]}
                                       for k in range(len(means))]
        return

    # fitting
    def _fit_cl_dna_experiment(self):
        self.coefficients['p_prot'] = self.parameters['initial_guess']['p_prot']
        self.coefficients['p_dna'] = self.parameters['initial_guess']['p_dna']
        self.coefficients['k_dna'] = self.parameters['initial_guess']['k_dna']
        self.coefficients['k_prot'] = None
        self.coefficients['l_dna'] = self.parameters['initial_guess']['l_dna']

        self.data['x_dna'] = np.array([invert_wlc(f, self.coefficients['p_dna'], self.coefficients['k_dna'])
                                       for f in self.data['F']])
        self.data['d_dna'] = self.coefficients['l_dna'] * self.data['x_dna']
        return

    def _fit_cl_dna_theory(self):
        self.coefficients['p_prot'] = self.parameters['initial_guess']['p_prot']
        self.coefficients['p_dna'] = self.parameters['initial_guess']['p_dna']
        self.coefficients['k_dna'] = self.parameters['initial_guess']['k_dna']
        self.coefficients['k_prot'] = self.parameters['initial_guess']['k_prot']
        self.coefficients['l_dna'] = self.parameters['initial_guess']['l_dna']

        self.data['x_dna'] = np.array([invert_wlc(f, self.coefficients['p_dna'], self.coefficients['k_dna'])
                                       for f in self.data['F']])
        self.data['d_dna'] = self.coefficients['l_dna'] * self.data['x_dna']
        return

    def _fit_cl_none_experiment(self):
        self.coefficients['p_prot'] = self.parameters['initial_guess']['p_prot']
        self.coefficients['p_dna'] = None
        self.coefficients['k_dna'] = None
        self.coefficients['k_prot'] = None
        self.coefficients['l_dna'] = None

        self.data['x_dna'] = np.zeros(len(self.data))
        self.data['d_dna'] = np.zeros(len(self.data))
        return

    def _fit_cl_none_theory(self):
        self.coefficients['p_prot'] = 0.7735670704545268
        self.coefficients['p_dna'] = self.parameters['initial_guess']['p_dna']
        self.coefficients['k_dna'] = self.parameters['initial_guess']['k_dna']
        self.coefficients['k_prot'] = 200
        self.coefficients['l_dna'] = self.parameters['initial_guess']['l_dna']

        self.data['x_dna'] = np.zeros(len(self.data))
        self.data['d_dna'] = np.zeros(len(self.data))
        return

    def fit_contour_lengths(self):
        if self.logger:
            self.logger.info("Fitting contour lengths...")
        if self.parameters['linker'] == 'dna' and self.parameters['source'] == 'experiment':
            self._fit_cl_dna_experiment()
        elif self.parameters['linker'] == 'dna' and self.parameters['source'] == 'theory':
            self._fit_cl_dna_theory()
        elif self.parameters['linker'] == 'none' and self.parameters['source'] == 'experiment':
            self._fit_cl_none_experiment()
        elif self.parameters['linker'] == 'none' and self.parameters['source'] == 'theory':
            self._fit_cl_none_theory()
        else:
            raise ValueError("Unknown combination of data source and linker. \n"
                             "Got data source " + self.parameters['source'] + " and linker " +
                             self.parameters['linker'])

        self.data['x_prot'] = np.array([invert_wlc(f, self.coefficients['p_prot'], self.coefficients['k_prot'])
                                        for f in self.data['F']])
        self.data['d_prot'] = self.data['d'] - self.data['d_dna']
        self.data['hist_values'] = self.data['d_prot'] / self.data['x_prot']
        self._parse_cl_histogram()
        return

    def find_rupture_forces(self):
        if self.logger:
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
        if self.logger:
            print("\t-> " + str(datetime.now().time()) + "\t Finding energies...")
        return

    # Plotting #
    def _plot_histogram(self, position=None, max_contour_length=None):
        position.set_title('Contour length histogram (trace ' + str(self.name) + ')')
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Counts')
        if max_contour_length:
            position.set_xlim(0, max_contour_length)
        lspace = np.linspace(self.ranges[0][0], self.ranges[-1][1])

        for state in range(len(self.ranges)):
            current_range = self.ranges[state]
            mu = self.coefficients['l_prot'][state]['mean']
            gamma = self.coefficients['l_prot'][state]['width']
            values = [v for v in self.data['hist_values'] if current_range[0] <= v <= current_range[1]]
            residues = 1 + int(mu / self.parameters['residues_distance'])
            label = "L= " + str(round(mu, 3)) + ' (' + str(residues) + ' AA)'
            n, bins, patches = position.hist(values, bins=200, color=colors[state % len(colors)], alpha=0.5)
            area = find_area(n, bins)
            position.plot(lspace, area * norm.pdf(lspace, mu, gamma), linestyle='--', label=label,
                          color=colors[state % len(colors)])
        position.legend()
        return

    def _plot_fits(self, position=None):
        position.set_title('Trace fits (trace ' + str(self.name) + ')')
        position.set_xlim(min(self.data['d']), max(self.data['d']))
        position.set_ylim(0, max(self.data['F']))
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Force [pN]')
        position.plot(self.data['d'], self.data['F'], label='Original data')
        position.plot(self.smoothed['d'], self.smoothed['F'], label='Smoothed data')

        p_prot = self.coefficients.get('p_prot', 0)
        p_dna = self.coefficients.get('p_dna', 0)
        k_dna = self.coefficients.get('k_dna', None)
        k_prot = self.coefficients.get('k_prot', None)
        l_dna = self.coefficients.get('l_dna', 0)

        f_plot = np.linspace(0.1, max(self.data['F']))
        if p_dna > 0:
            d_dna = l_dna * np.array([invert_wlc(f, p_dna, k_dna) for f in f_plot])
        else:
            d_dna = np.zeros(len(f_plot))
        xprot = np.array([invert_wlc(f, p_prot, k_prot) for f in f_plot])

        for state in range(len(self.coefficients['lengths'])):
            lprot = self.coefficients['l_prot'][state]['mean']
            residues = 1 + int(lprot / self.parameters['residues_distance'])
            d_prot = lprot * xprot
            d_plot = d_dna + d_prot
            label = 'Fit L=' + str(round(lprot, 3)) + ' (' + str(residues) + ' AA)'
            position.plot(d_plot, f_plot, label=label, color=colors[state % len(colors)])
            x0, x1 = self.ranges[state]
            position.hlines(max(self.data['F']) / 2, x0 + 2, x1 - 2, color=colors[state % len(colors)], alpha=0.5)
            position.axvline(x=x0, linestyle='--', color='#808080', linewidth=0.5)
            position.axvline(x=x1, linestyle='--', color='#808080', linewidth=0.5)

        position.legend()
        return

    # analyzing
    def analyze(self):
        if self.logger:
            self.logger.info("Analyzing trace " + str(self.name))
        self._generate_smooths()
        self.fit_contour_lengths()
        # self.find_rupture_forces()
        # self.find_energies()
        if self.logger:
            print(str(datetime.now().time()) + "\t Finished trace " + str(self.name))
        return

    def summary(self):
        return

    def get_info(self):
        info = []
        info.append("Trace name " + str(self.name))
        info.append("Trace source file " + str(self.orig_input))
        return '\n'.join(info)
