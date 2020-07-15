from .stretching_tools import *
from .default_parameters import default_parameters
from .trace import Trace
from matplotlib import pyplot as plt
import os


class Structure:
    """The class containing information about the whole experiment, consisting of many measurements (traces).

    Attributes:
        input_data (str/Pandas Dataframe/None): The data to be analyzed or the path the file with data.
        parameters (dict): The dictionary of measurement parameters.
        rupture_forces (dict): The dictionary with measured rupture forces.

    """
    def __init__(self, input_data=None, cases=None, columns=None, name=None, debug=False, **kwargs):
        self.orig_input = input_data
        self.debug = debug

        # setting the name
        self._set_name(input_data, name)

        # initializing log file if debug
        if self.debug:
            logger = set_logger(self.name, mode='w')
            logger.debug('Initializing the "Structure" class with name: ' + str(self.name))
            close_logs(logger)

        # setting parameters
        self.parameters = {}
        self._set_parameters(**kwargs)

        # reading the data (setting the traces)
        self.traces = []
        self._read_data(input_data, cases, columns)
        return

    ''' Methods for users '''
    def add_trace(self, filename, cases=None, columns=None, **kwargs):
        """ Adding the trace to the structure.

                Args:
                    filename (str/Pandas Dataframe): Path to the file, or Pandas Dataframe with the data of the trace.
                    cases (list, optional): The list of cases to be added. If None, all cases will traces will be used.
                        Default: None.
                    columns ([str, str], optional): The list of headers of the columns containing the distance and
                        corresponding force. If None, first two columns will be used. Default: None
                    **kwargs: The parameters which are going to be set for the traces. Accepted parameters are:
                    'linker' ('dna'/None), 'source' ('experiment'/'cg'/'aa'), 'speed' (float), 'states' (int),
                    'residues_distance', (float), 'low_force_cutoff' (float), 'initial_guess' (dictionary), 'sheet_name'
                    (string), 'separator' (char), unit (string).

                Returns:
                    List of the names of the Traces added.
                """

        n = len(self.traces)
        self._read_data(filename, cases=cases, columns=columns, **kwargs)
        names = [str(x) for x in range(n, len(self.traces))]
        return names

    def set(self, **kwargs):
        """ Setting the parameter/coefficient to the Structure. The parameters to be set are given as keyword arguments.

                Args:
                    **kwargs: The parameters which are going to be set. Accepted parameters are: 'linker' ('dna'/None),
                    'source' ('experiment'/'cg'/'aa'), 'speed' (float), 'states' (int), 'residues_distance', (float),
                    'low_force_cutoff' (float), 'initial_guess' (dictionary), 'plot_columns' (int).

                Returns:
                    True if successful, False otherwise.
                """

        self.parameters = {**self.parameters, **kwargs}
        for trace in self.traces:
            trace.set(**kwargs)

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Structure parameters updated. Current parameters: " + str(self.parameters))
            close_logs(logger)
        return True

    def make_plots(self):
        """ Plotting all the output figures.


                Returns:
                    True if successful, False otherwise.
                """

        self.make_histograms()
        self.make_partial_plots()
        return True

    def make_partial_plots(self, output=None):
        """ Plot the fits for each trace in structure.

                Args:
                    output (str, optional): The output name for the figure.

                Returns:
                    True if successful, False otherwise.
                """

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Plotting fits of individual traces")
            close_logs(logger)

        # setting the figure
        number = len(self.traces)                                   # number of traces to plot
        columns = self.parameters['plot_columns']                   # number of columns
        rows = max(int(np.ceil(float(number) / columns)), 2)        # number of rows
        fig = plt.subplots(dpi=600, figsize=(10*int(columns), 5*int(rows)))
        axes = []
        max_contour_length = self.parameters['boundaries'][1]     # to align plots

        # plotting each trace
        for k in range(0, number):
            axes.append(plt.subplot2grid((rows, 2 * columns), (int(k / columns), (2 * k) % (2 * columns))))
            self.traces[k].plot_histogram(position=axes[-1], max_contour_length=max_contour_length)
            axes.append(plt.subplot2grid((rows, 2 * columns), (int(k / columns), ((2 * k) + 1) % (2 * columns))))
            self.traces[k].plot_fits(position=axes[-1])
        plt.tight_layout()

        if not output:
            ofile = self.name + '_contour_lengths.png'
        else:
            ofile = output

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Saving fits figure to " + ofile)
            close_logs(logger)

        plt.savefig(ofile)
        plt.close(fig)
        return True

    def make_histograms(self, output=None):
        """ Plot the cumulative histograms for all traces in structure.

                Args:
                    output (str, optional): The output name for the figure.

                Returns:
                    True if successful, False otherwise.
                """

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Making histograms.")
            close_logs(logger)

        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        self._plot_total_contour_length_histo(axes[0, 0])   # the total contour length histogram
        self._plot_overlaid_traces(axes[0, 1])              # the overlaid smoothed traces
        self._plot_forces_histogram(axes[1, 0])             # the rupture forces histogram
        self._plot_dhs_analysis(axes[1, 1])                 # Dudko analysis
        fig.tight_layout()

        if not output:
            ofile = self.name + '_histograms.png'
        else:
            ofile = output

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Saving histograms figure to " + ofile)
            close_logs(logger)

        plt.savefig(ofile)
        plt.close(fig)
        return True

    def plot_contour_length(self, output=None):
        """ Plot histogram of contour length.

                 Args:
                     output (str, optional): The output name for the figure.

                 Returns:
                     True if successful, False otherwise.
                 """

        fig, axes = plt.subplots(1, 1, dpi=600, figsize=(5, 5))
        self._plot_total_contour_length_histo(axes[0])
        if not output:
            ofile = self.name + '_contour_length.png'
        else:
            ofile = output

        plt.savefig(ofile)
        plt.close()
        return

    def plot_traces(self, output=None):
        """ Plot overlaid traces.

                  Args:
                      output (str, optional): The output name for the figure.

                  Returns:
                      True if successful, False otherwise.
                  """

        fig, axes = plt.subplots(1, 1, dpi=600, figsize=(5, 5))
        self._plot_overlaid_traces(axes[0])

        if not output:
            ofile = self.name + '_traces.png'
        else:
            ofile = output

        plt.savefig(ofile)
        plt.close()
        return True

    def plot_forces(self, output=None):
        """ Plot histogram of rupture forces.

                  Args:
                      output (str, optional): The output name for the figure.

                  Returns:
                      True if successful, False otherwise.
                  """

        fig, axes = plt.subplots(1, 1, dpi=600, figsize=(5, 5))
        self._plot_forces_histogram(axes[0])

        if not output:
            ofile = self.name + '_rupture_forces.png'
        else:
            ofile = output

        plt.savefig(ofile)
        plt.close()
        return True

    def plot_dhs(self, output=None):
        """ Plot Dudko-Hummer-Szabo analysis.

                  Args:
                      output (str, optional): The output name for the figure.

                  Returns:
                      True if successful, False otherwise.
                  """

        fig, axes = plt.subplots(1, 1, dpi=600, figsize=(5, 5))
        self._plot_dhs_analysis(axes[0])

        if not output:
            ofile = self.name + '_dhs_analysis.png'
        else:
            ofile = output

        plt.savefig(ofile)
        plt.close()
        return True

    def collect_coefficients(self):
        """ Collecting the coefficients for the whole structure as the means of coefficients from all traces.

                  Returns:
                      True if successful, False otherwise.
                  """

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Collecting coefficients")
            close_logs(logger)

        # collecting the traces coefficients
        coefficients = ['p_prot', 'k_prot', 'p_dna', 'k_dna', 'l_dna']
        for c in coefficients:
            self.parameters[c] = np.mean(np.array([t.parameters.get(c, 0) for t in self.traces]))

        # creating the data for contour length histogram
        for t in self.traces:
            self.hist_values += list(t.data['hist_values'])

        # decomposing the contour length histogram = obtaining the contour lengths
        parameters, boundaries = decompose_histogram(np.array(self.hist_values), states=self.parameters['states'],
                                                     significance=self.parameters['significance'])
        self.parameters['l_prot'] = parameters
        self.parameters['boundaries'] = boundaries

        # collecting the rupture forces
        self._analyze_rupture_forces()

        # perform Dudko-Hummer-Szabo analysis
        self._analyze_dhs()
        return True

    def analyze(self):
        """ Perform the whole analysis.

                  Returns:
                      True if successful, False otherwise.
                  """
        if self.debug:
            logger = set_logger(self.name)
            logger.info("Starting analysis of the experiment")
            close_logs(logger)

        for t in self.traces:
            t.analyze()

        self.collect_coefficients()
        self.make_plots()
        self.save_data()
        return

    def save_data(self, output=None):
        """ Save the fitted data to the text file.

                  Args:
                      output (str, optional): The output name for the file.

                  Returns:
                      True if successful, False otherwise.
                  """

        separator = '################\n'
        if not output:
            oname = str(self.name) + '_results'
        else:
            oname = output

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Saving data to " + oname)
            close_logs(logger)

        result = [self._get_general_info(separator),
                  self._get_traces_info(separator),
                  self._get_cummulative_statistics(separator),
                  self._get_force_analysis_info(separator)]

        with open(oname, 'w') as ofile:
            ofile.write('\n'.join(result))
        return True

    ''' Protected methods '''
    # getting info
    def _get_general_info(self, separator):
        info = ['Experiment name\t' + str(self.name),
                'Experiment source file\t' + str(self.orig_input),
                'Number of traces:\t' + str(len(self.traces)),
                'Data source\t' + str(self.parameters['source']),
                'Data unit\t' + str(self.parameters['unit']),
                'Structure linker\t' + str(self.parameters['linker']),
                'Number of states\t' + str(self.parameters['states']),
                'Pulling speed\t' + str(self.parameters['speed']),
                'Low force cutoff\t' + str(self.parameters['low_force_cutoff']),
                'Significance\t' + str(self.parameters['significance']),
                'Fit initial guess\t' + str(self.parameters['initial_guess']),
                'Residue distance (nm)\t' + str(self.parameters['residue_distance']),
                'Sheet name\t' + str(self.parameters['sheet_name']),
                'Separator\t' + str(self.parameters['separator']),
                separator]
        return '\n'.join(info)

    def _get_traces_info(self, separator):
        result = ['Summary of individual curves'] + [t.get_info() for t in self.traces] + [separator]
        return '\n'.join(result)

    def _get_cummulative_statistics(self, separator):
        result = ['Summary of the cumulative statistics',
                  'p_Prot:\t\t' + str(self.parameters['p_prot']),
                  'k_Prot:\t\t' + str(self.parameters['k_prot']),
                  'p_DNA:\t\t' + str(self.parameters['p_dna']),
                  'k_DNA:\t\t' + str(self.parameters['k_dna']),
                  'l_DNA:\t\t' + str(self.parameters['l_dna']),
                  'Contour length\tgamma\t',
                  self.parameters['l_prot'].to_csv(sep='\t'),
                  separator]
        return '\n'.join(result)

    def _get_force_analysis_info(self, separator):
        result = ['Dudko-Hummer-Szabo analysis', self.forces.to_csv(sep='\t')]
        result += [str(key) + '\t' + str(v) + '\t' + str(self.dhs_results[key][v])
                   for key in self.dhs_results.keys() for v in self.dhs_results[key].keys()]
        result.append(separator)
        return '\n'.join(separator)

    # reading
    def _read_data(self, input_data=None, cases=None, columns=None, **kwargs):
        """ Reading the data and splitting them into Traces """
        parameters = {**self.parameters, **kwargs}

        if isinstance(input_data, pd.DataFrame):
            data = read_dataframe(input_data, cases, columns)
        else:
            data = read_from_file(input_data, cases, columns, parameters, self.name, self.debug)

        # preprocessing data
        data.dropna(axis='columns', how='all', inplace=True)
        if len(data) == 0:
            if self.debug:
                logger = set_logger(self.name)
                logger.info("Initializing empty class. Hope you'll add some traces to analyze.")
                close_logs(logger)
                return

        # dividing data into traces
        headers = list(data)
        for k in range(len(headers) % 2, len(headers), 2):    # if the first column is the index
            if headers[k][0] == 'd':
                trace_name = headers[k].strip('d_')
            else:
                trace_name = str(int(k/2))
            trace_data = data[[headers[k], headers[k + 1]]]
            new_headers = {headers[k]: 'd', headers[k + 1]: 'F'}
            trace_data = trace_data.rename(columns=new_headers)
            if parameters['unit'] == 'A':
                trace_data['d'] = trace_data['d'] / 10
            self.traces.append(Trace(trace_name, trace_data, input_data, experiment_name=self.name, debug=self.debug,
                                     parameters=self.parameters))
        return True

    # setting
    def _set_name(self, filename, name):
        """ Setting the Structure name """
        if name:
            self.name = name
        elif isinstance(filename, str) and len(filename) > 0:       # possibly the path to file
            self.name = os.path.splitext(os.path.basename(filename))[0]
        else:
            self.name = 'Experiment'
        return 0

    def _set_parameters(self, **kwargs):
        """ The method of setting the parameters for the whole Structure """
        # filling in the default parameters
        self.parameters = default_parameters

        # substituting parameters
        for key in kwargs:
            self.parameters[key] = kwargs[key]

        # dealing with initial_guess
        if 'initial_guess' not in list(kwargs.keys()):
            self.parameters['initial_guess'] = self.parameters['initial_guess'][self.parameters['source']]

        # setting other attributes
        self.hist_values = []
        self.forces = pd.DataFrame(columns=['means', 'rupture_forces'])
        self.states = []
        self.max_f = 0
        self.dhs_states = []
        self.dhs_results = {}

        if self.debug:
            logger = set_logger(self.name)
            logger.debug("Set structure parameters:\t" + str(self.parameters))
            close_logs(logger)
        return

    # plotting
    """ Plotting the contour length panel """
    def _plot_total_contour_length_histo(self, position):
        if not self.hist_values:
            if self.debug:
                logger = set_logger(self.name)
                logger.debug("Structure.hist_values empty, nothing to plot. Did you run .analyze() or "
                             ".collect_coefficients ?")
                close_logs(logger)
            return
        if self.parameters['boundaries']:
            bound = self.parameters['boundaries'][1]
        else:
            bound = max(self.hist_values)

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Making cumulative contour length histograms.")
            close_logs(logger)

        # setting the scene
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')
        position.set_xlim(0, bound)

        # plotting histogram
        position.hist(self.hist_values, bins=500, density=True, alpha=0.5)

        # plotting decomposition
        plot_decomposed_histogram(position, self.parameters['l_prot'], bound, self.parameters['residues_distance'])

        position.legend()
        return

    def _plot_overlaid_traces(self, position):
        if self.debug:
            logger = set_logger(self.name)
            logger.info("Plotting overlaid traces.")
            close_logs(logger)

        # setting the scene
        position.set_xlabel('Distension [nm]')
        position.set_ylabel('Force [pN]')
        position.set_title('All smoothed traces overlaid')
        max_f = 0
        position.set_ylim(0, max_f)

        # plotting overlaid smoothed traces
        for k in range(len(self.traces)):
            t = self.traces[k]
            position.plot(t.smoothed['d'], t.smoothed['F'], color=get_gray_shade(k, len(self.traces)))
            max_f = max(max_f, t.data['F'].max())

        # plotting fits
        plot_trace_fits(position, self.parameters, max_f, self.parameters['residues_distance'])

        position.legend()
        return

    def _plot_forces_histogram(self, position):
        if not self.forces:
            if self.debug:
                logger = set_logger(self.name)
                logger.debug("Structure.forces empty, nothing to plot. Did you run .analyze() or "
                             ".collect_coefficients ?")
                close_logs(logger)
            return

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Making the rupture forces histogram.")
            close_logs(logger)

        # setting the scene
        position.set_xlabel('Force [pN]')
        position.set_ylabel('Occurences')
        position.set_title('Rupture force histogram')
        f_space = np.linspace(0, self.max_f)

        # plotting the histograms and the fitted contour
        for k in range(len(self.states)):
            state = self.states[k]
            data_to_plot = np.array(self.forces[self.forces['state'] == state]['rupture_forces'].dropna())
            # TODO the line below is redundant with _analyze_dhs
            parameters, boundaries = decompose_histogram(data_to_plot, significance=self.parameters['significance'])

            y_plot = np.zeros(len(f_space))
            for index, row in parameters.iterrows():    # if there is more than 1 contour length for this force
                y_plot += gauss(f_space, row['heights'], row['means'], row['widths'])

            label = '; '.join([str(round(x, 3)) for x in list(parameters['means'].values)]) + ' pN'
            position.hist(data_to_plot, bins=20, color=get_color(k, len(self.states)), alpha=0.5, density=True,
                          label=label)
            position.plot(f_space, y_plot, linestyle='--', color=get_color(k, len(self.states)))

        position.legend()
        return

    def _plot_dhs_analysis(self, position):
        # setting the scene
        position.set_title('Dudko-Hummer-Szabo lifetime')
        position.set_xlabel('Rupture force [pN]')
        position.set_ylabel('State lifetime [s]')
        position.set_ylim(top=1000)
        position.set_yscale('log')

        for k in range(len(self.dhs_states)):
            beg, end, l_prot, force, lifetime = self.dhs_states[k]
            f_space = np.linspace(beg, end, 100)
            label = 'l_prot=' + str(round(l_prot, 3)) + '; force=' + str(round(force, 3))
            position.plot(f_space, lifetime, label=label, color=get_color(k, len(self.dhs_states)))
        position.legend()
        return

    # analyzing
    def _analyze_rupture_forces(self):
        """ Collecting and analyzing rupture forces"""

        # collecting forces
        self.forces = pd.concat([t.parameters['l_prot'][['means', 'rupture_forces']].dropna() for t in self.traces],
                                ignore_index=True)
        self.max_f = self.forces['rupture_forces'].max()

        # assigning a rupture force to the state from the list of cumulative states
        # TODO - calculate it using Cauchy distribution
        for index, row in self.parameters['l_prot'][['means', 'widths', 'heights']].iterrows():
            mean, width, height = tuple(row.to_numpy())
            self.forces['state_' + str(index)] = gauss(np.array(self.forces['means']), height, mean, width)
            self.states.append('state_' + str(index))

        # selecting the corresponding state
        self.forces['state'] = self.forces[self.states].idxmax(axis=1)
        return True

    def _analyze_dhs(self):
        # TODO clean it up
        """ Calculating the data for the Dudko-Hummer-Szabo analysis of the states lifetime"""
        for k in range(len(self.states)):
            state = self.states[k]
            data = np.array(self.forces[self.forces['state'] == state]['rupture_forces'].dropna())
            parameters, boundaries = decompose_histogram(data, significance=self.parameters['significance'])
            state_index = int(state.strip('state_'))
            l_prot = self.parameters['l_prot'].loc[state_index, 'means']
            beg, end = boundaries
            f_space = np.linspace(beg, end, 100)

            for ind, row in parameters[['means', 'widths', 'heights']].iterrows():
                mean, width, height = tuple(row.to_numpy())

                force_load = get_force_load(f_space, self.parameters['speed'], l_prot,
                                            self.parameters['p_prot'], self.parameters['k_prot'])
                probability = gauss(f_space, height, mean, width)
                denominator = probability * force_load
                nominator = integrate_gauss(f_space, mean, width)

                # calculating the lifetime and fitting it
                lifetime = nominator / denominator
                self.dhs_states.append([beg, end, l_prot, mean, lifetime])
                self.dhs_results[(round(l_prot, 3), round(mean, 3))] = dhs_feat(f_space, lifetime)
        return True
