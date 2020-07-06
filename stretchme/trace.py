from datetime import datetime
from matplotlib import pyplot as plt
from .stretching_tools import *
import pandas as pd


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
        if self.logger:
            self.logger.info("Trace name " + str(name) + " got and parsed the data.")

        self.smoothed = pd.DataFrame()

        self.boundaries = []
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

    def _parse_cl_histogram(self, significance=0.001):
        parameters, boundaries = decompose_histogram(np.array(self.data['hist_values']), significance)
        self.boundaries = boundaries
        self.coefficients['l_prot'] = parameters
        return

    # setting
    def set_residues(self, value):
        self.parameters['residues'] = value
        if self.logger:
            self.logger.info("Set trace 'residues' value to " + str(value))
        return

    def set_distance(self, value):
        self.parameters['distance'] = value
        if self.logger:
            self.logger.info("Set trace 'distance' value to " + str(value))
        return

    def set_linker(self, value):
        self.parameters['linker'] = value
        if self.logger:
            self.logger.info("Set trace 'linker' value to " + str(value))
        return

    def set_source(self, value):
        self.parameters['source'] = value
        if self.logger:
            self.logger.info("Set trace 'source' value to " + str(value))
        return

    def set_speed(self, value):
        self.parameters['speed'] = value
        if self.logger:
            self.logger.info("Set trace 'speed' value to " + str(value))
        return

    def set_residues_distance(self, value):
        self.parameters['residues_distance'] = value
        if self.logger:
            self.logger.info("Set trace 'residues_distance' value to " + str(value))
        return

    def set_minimal_stretch_distance(self, value):
        self.parameters['minimal_stretch_distance'] = value
        if self.logger:
            self.logger.info("Set trace 'minimal_stretch_distance' value to " + str(value))
        return

    def set_high_force_cutoff(self, value):
        self.parameters['high_force_cutoff'] = value
        if self.logger:
            self.logger.info("Set trace 'high_force_cutoff' value to " + str(value))
        return

    def set_low_force_cutoff(self, value):
        self.parameters['low_force_cutoff'] = value
        if self.logger:
            self.logger.info("Set trace 'low_force_cutoff' value to " + str(value))
        return

    def set_max_rupture_force(self, value):
        self.parameters['max_rupture_force'] = value
        if self.logger:
            self.logger.info("Set trace 'max_rupture_force' value to " + str(value))
        return

    def set_initial_guess(self, value):
        self.parameters['initial_guess'] = value
        if self.logger:
            self.logger.info("Set trace 'initial_guess' value to " + str(value))
        return

    def set_states(self, value):
        self.parameters['states'] = value
        if self.logger:
            self.logger.info("Set trace 'states' value to " + str(value))
        return

    # fitting
    def _fit_distances(self):
        for index, row in self.coefficients['l_prot'].iterrows():
            l_prot = row['means']
            self.smoothed['state_' + str(index)] = np.array([wlc(d, l_prot, self.coefficients.get('p_prot', 0),
                                            self.coefficients.get('k_prot', None)) for d in list(self.smoothed['d'])])
        return

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
        # TODO finish it
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
        # TODO finish it
        self.coefficients['p_prot'] = self.parameters['initial_guess']['p_prot']
        self.coefficients['p_dna'] = 0
        self.coefficients['k_dna'] = None
        self.coefficients['k_prot'] = None
        self.coefficients['l_dna'] = 0

        self.data['x_dna'] = np.zeros(len(self.data))
        self.data['d_dna'] = np.zeros(len(self.data))
        return

    def _fit_cl_none_theory(self):
        self.coefficients['k_dna'] = None
        self.coefficients['l_dna'] = 0
        self.coefficients['p_dna'] = 0
        self.data['x_dna'] = np.zeros(len(self.data))
        self.data['d_dna'] = np.zeros(len(self.data))


        p_prot = self.parameters['initial_guess']['p_prot']
        k_prot = self.parameters['initial_guess']['k_prot']
        p_prot, k_prot = minimize_pk(self.data[['d', 'F']], self.smoothed, p_prot, k_prot)
        print(p_prot, k_prot)
        self.coefficients['p_prot'] = p_prot
        self.coefficients['k_prot'] = k_prot
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
        self._fit_distances()
        return

    def find_rupture_forces(self, max_distance=0.3):
        # the ranges are determined by the proximity of fitted curve to the smoothed data
        # TODO clean it up
        if self.logger:
            self.logger.info("Finding rupture forces...")
        last_end = 0
        forces = []
        for ind, row in self.coefficients['l_prot'].iterrows():
            data_close = self.smoothed[abs(self.smoothed['F'] - self.smoothed['state_' + str(ind)]) <= max_distance]
            data_close = data_close[data_close['d'] > last_end]['d']
            beg = data_close.min()
            end = data_close.max()
            last_end = end
            forces.append(self.data[self.data['d'].between(beg, end)]['F'].max())

        self.coefficients['l_prot']['rupture_forces'] = np.array(forces)
        return

    def find_energies(self):
        if self.logger:
            self.logger.info("Finding energies...")
        return

    # Plotting #
    def _plot_histogram(self, position=None, max_contour_length=None):
        position.set_title('Contour length histogram (trace ' + str(self.name) + ')')
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Counts')
        if not max_contour_length:
            if self.boundaries:
                max_contour_length = self.boundaries[1]
            else:
                max_contour_length = self.data['hist_values'].max()
        position.set_xlim(0, max_contour_length)
        l_space = np.linspace(0, max_contour_length)

        # whole histogram
        position.hist(self.data['hist_values'], bins=500, density=True, alpha=0.5)

        # decomposed histogram
        plot_decomposed_histogram(position, self.coefficients['l_prot'], l_space, self.parameters['residues_distance'])

        position.legend()
        return

    def _plot_fits(self, position=None):
        position.set_title('Trace fits (trace ' + str(self.name) + ')')
        position.set_xlim(min(self.data['d']), max(self.data['d']))
        position.set_ylim(0, self.data['F'].max())
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Force [pN]')

        #plotting original data
        position.plot(self.data['d'], self.data['F'], label='Original data')
        position.plot(self.smoothed['d'], self.smoothed['F'], label='Smoothed data')

        # plotting fits
        f_space = np.linspace(0.1, self.data['F'].max())
        plot_trace_fits(position, self.coefficients, f_space, self.parameters['residues_distance'])

        position.legend()
        return

    # analyzing
    def analyze(self):
        print("Analyzing " + str(self.name))
        if self.logger:
            self.logger.info("Analyzing trace " + str(self.name))
        self._generate_smooths()
        self.fit_contour_lengths()
        self.find_rupture_forces()
        # self.find_energies()
        if self.logger:
            self.logger.info("Finished trace " + str(self.name))
        return

    def summary(self):
        result = []
        separator = '####\n'
        result.append(str(self.name))
        result.append('->\tp_Prot:\t\t' + str(self.coefficients['p_prot']))
        result.append('->\tp_DNA:\t\t' + str(self.coefficients['p_dna']))
        result.append('->\tk_Prot:\t\t' + str(self.coefficients['k_prot']))
        result.append('->\tk_DNA:\t\t' + str(self.coefficients['k_dna']))
        result.append('->\tl_DNA:\t\t' + str(self.coefficients['l_dna']))
        result.append('->\tContour length\tgamma\tsigma^2\trupture_force')
        result.append('->\t' + self.coefficients['l_prot'].to_csv(sep='\t'))
        result.append('->\tContour length boundaries')
        result.append('->\t' + str(self.boundaries))
        result.append(separator)

        return '\n'.join(result)

    def get_info(self):
        info = []
        info.append("Trace name " + str(self.name))
        info.append("Trace source file " + str(self.orig_input))
        return '\n'.join(info)
