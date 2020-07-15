from .stretching_tools import *
import pandas as pd
from matplotlib import pyplot as plt


class Trace:
    """The class containing information about one, single measurement.

    Attributes:
        name (str): The class name.
        orig_input (str): The original input (e.g. the path to the source file).
        parameters (dict): The dictionary of measurement parameters and fitted coefficients.
        data (Pandas Dataframe): The dataframe with columns 'd' and 'F' corresponding to measured distance and force.
        rupture_forces (dict): The dictionary with measured rupture forces.

    """
    def __init__(self, name, data, orig_input, experiment_name='Experiment', debug=False, parameters=None):
        self.name = str(name)
        self.orig_input = orig_input
        self.experiment_name = experiment_name
        self.debug = debug

        # setting the trace parameters
        if not parameters:
            self.parameters = default_parameters
            self.parameters['initial_guess'] = default_parameters['initial_guess'][None]
        else:
            self.parameters = parameters

        # filling the data
        self.data = data.dropna()
        self.data = self.data[self.data['F'] > self.parameters['low_force_cutoff']].sort_values(by='d')

        self.smoothed = pd.DataFrame()
        self.boundaries = []
        self.state_boundaries = pd.DataFrame(columns=['state', 'beg', 'end'])
        self.rupture_forces = {}

        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Trace name " + str(name) + " got and parsed the data.")
            close_logs(logger)
        return

    ''' Methods for users '''
    def set(self, **kwargs):
        """ Method of setting the parameter to the trace.

                Args:
                    **kwargs: The parameters which are going to be set. Accepted parameters are: 'linker' ('dna'/None),
                    'source' ('experiment'/'cg'/'aa'), 'speed' (float), 'states' (int), 'residues_distance', (float),
                    'low_force_cutoff' (float), and 'initial_guess' (dictionary).

                Returns:
                    True if successful, False otherwise.
                """
        self.parameters = {**self.parameters, **kwargs}
        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Trace parameters updated. Current parameters: " + str(self.parameters))
            close_logs(logger)
        return True

    def analyze(self):
        """ Performing the whole analysis of the trace.

                Returns:
                    True if successful, False otherwise.
                """
        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Starting analysis of the trace.")
            close_logs(logger)

        self.generate_smooths()
        self.fit_contour_lengths()
        self.find_ranges()
        self.find_rupture_forces()

        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Finished trace " + str(self.name))
            close_logs(logger)

        return True

    def fit_contour_lengths(self):
        """ Fitting the contour lengths for the trace.

                Returns:
                    True if successful, False otherwise.
                """

        # fitting coefficients
        coefficients = fit_coefficients(self.data, self.parameters)
        self.parameters = {**self.parameters, **coefficients}

        # transforming coordinates
        self.data['x_prot'] = np.array([invert_wlc(f, self.parameters['p_prot'], self.parameters['k_prot'])
                                        for f in self.data['F']])
        self.data['d_dna'] = get_d_dna(self.parameters['p_dna'], self.parameters['l_dna'], self.parameters['k_dna'],
                                       self.data['F'])
        self.data['d_prot'] = self.data['d'] - self.data['d_dna']
        self.data['hist_values'] = self.data['d_prot'] / self.data['x_prot']

        # decomposing contour length histogram
        parameters, bounds = decompose_histogram(np.array(self.data['hist_values']), states=self.parameters['states'],
                                                 significance=self.parameters['significance'])
        self.boundaries = bounds
        self.parameters['l_prot'] = parameters
        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Contour length fitted. Coefficient got:\n" + str(coefficients))
            close_logs(logger)
        return True

    def find_ranges(self, max_distance=None):
        """ Method of setting the parameter to the trace.

                Args:
                   max_distance (float/None, optional): The maximal deviation of the fitted curve from the data in the
                    state. If None, the default value will be used.

                Returns:
                    True if successful, False otherwise.
                """
        if not max_distance:
            max_distance = self.parameters['max_distance']
        last_end = 0
        for index, row in self.parameters['l_prot'].iterrows():
            state = 'state_' + str(index)
            l_prot = row['means']
            self.smoothed[state] = \
                np.array([ewlc(d, l_prot, self.parameters.get('p_prot', 0), self.parameters.get('k_prot', None))
                          for d in list(self.smoothed['d'])])
            data_close = self.smoothed[abs(self.smoothed['F'] - self.smoothed[state]) <= max_distance]
            data_close = data_close[data_close['d'] > last_end]['d']
            last_end = data_close.max()
            new_row = pd.DataFrame({'state': state, 'beg': data_close.min(), 'end': last_end}, index=[0])
            self.state_boundaries = pd.concat([self.state_boundaries, new_row], ignore_index=True)

        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Ranges found. These are:\n" + str(self.state_boundaries))
            close_logs(logger)
        return True

    def find_rupture_forces(self):
        """ Find the rupture forces for each state.

                Returns:
                    True if successful, False otherwise.
                """

        forces = []
        for ind, row in self.state_boundaries.iterrows():
            state, beg, end = tuple(row.to_numpy())
            forces.append(self.data[self.data['d'].between(beg, end)]['F'].max())
        self.parameters['l_prot']['rupture_forces'] = np.array(forces)

        if self.debug:
            logger = set_logger(self.experiment_name)
            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_rows', None)
            logger.info("Rupture forces found. These are:\n" + str(self.parameters['l_prot']))
            pd.reset_option('display.max_columns')
            pd.reset_option('display.max_rows')
            close_logs(logger)
        return True

    def plot(self, output=None):
        """ Plots the contour length histogram and the fits.

                Args:
                    output (string/None, optional): The name of the output file.

                Returns:
                    True if successful, False otherwise.
                """
        fig, axes = plt.subplots(1, 2, dpi=600, figsize=(10, 5))
        self.plot_histogram(position=axes[0])
        self.plot_fits(position=axes[1])
        plt.tight_layout()

        if not output:
            ofile = self.experiment_name + '_' + self.name + '.png'
        else:
            ofile = output

        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Saving figure to: " + ofile)
            close_logs(logger)

        plt.savefig(ofile)
        plt.close(fig)
        return

    def get_info(self):
        """ Get the information about the trace.

                Returns:
                    String: The formatted information
                """
        separator = '####'
        info = ['Trace name ' + str(self.name),
                'Trace source file ' + str(self.orig_input),
                'p_Prot:\t\t' + str(self.parameters['p_prot']),
                'k_Prot:\t\t' + str(self.parameters['k_prot']),
                'p_DNA:\t\t' + str(self.parameters['p_dna']),
                'k_DNA:\t\t' + str(self.parameters['k_dna']),
                'l_DNA:\t\t' + str(self.parameters['l_dna']),
                'Contour length\tgamma\t',
                self.parameters['l_prot'].to_csv(sep='\t'),
                separator]
        return '\n'.join(info)

    def generate_smooths(self, intervals=None):
        """ Generate smooth trace, stored as trace.smoothed dataframe.

                Args:
                    intervals (int/None, optional): The number of distance intervals in the smoothed trace. If None, the
                        default value will be used.

                Returns:
                    True if successful, False otherwise.
                """
        if self.debug:
            logger = set_logger(self.experiment_name)
            logger.info("Generating smoothed traces.")
            close_logs(logger)

        if not intervals:
            intervals = self.parameters['intervals']

        f_smooth = []
        d_smooth = []

        # averaging over the interval
        d_range = np.linspace(min(self.data['d']), max(self.data['d']), intervals)
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
        return True

    ''' Restricted methods '''
    def plot_histogram(self, position=None, max_contour_length=None):
        """ Plotting the contour length histogram in a given position"""

        # setting the scene
        position.set_title('Contour length histogram (trace ' + str(self.name) + ')')
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Counts')
        if not max_contour_length:
            if self.boundaries:
                max_contour_length = self.boundaries[1]
            else:
                max_contour_length = self.data['hist_values'].max()
        position.set_xlim(0, max_contour_length)

        # whole histogram
        position.hist(self.data['hist_values'], bins=500, density=True, alpha=0.5)

        # decomposed histogram
        plot_decomposed_histogram(position, self.parameters['l_prot'], max_contour_length,
                                  self.parameters['residues_distance'])

        position.legend()
        return

    def plot_fits(self, position=None):
        """ Plotting the fits to the traces. """

        # setting the scene
        position.set_title('Trace fits (trace ' + str(self.name) + ')')
        position.set_xlim(min(self.data['d']), max(self.data['d']))
        position.set_ylim(0, self.data['F'].max())
        position.set_xlabel('Extension [nm]')
        position.set_ylabel('Force [pN]')

        # plotting original and smoothed data
        position.plot(self.data['d'], self.data['F'], label='Original data')
        position.plot(self.smoothed['d'], self.smoothed['F'], label='Smoothed data')

        # plotting fits
        plot_trace_fits(position, self.parameters, self.data['F'].max(), self.parameters['residues_distance'])

        position.legend()
        return
