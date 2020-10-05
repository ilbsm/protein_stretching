from .stretching_tools import *
from .default_parameters import default_parameters
from .trace import Trace
from matplotlib import pyplot as plt
import os
from copy import deepcopy


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
        print(self.parameters.get('state', '.') + '/' + self.name + '_' + str(n+1) + '_data.csv')
        if os.path.isfile(self.parameters.get('state', '.') + '/' + self.name + '_' + str(n+1) + '_data.csv'):
            trace_data = pd.read_csv(self.parameters.get('state', '.') + '/' + self.name + '_' + str(n+1) + '_data.csv')
            parameters = {**self.parameters, **kwargs}
            self.traces.append(Trace(self.name + '_' + str(n), trace_data, filename, experiment_name=self.name,
                                     debug=self.debug, parameters=parameters))
            print("Found and parsed data")
        else:
            self._read_data(filename, cases=cases, columns=columns, trace_name=self.name + '_' + str(n), **kwargs)
        names = [str(x) for x in range(n, len(self.traces))]
        return names

    def set(self, overide=True, **kwargs):
        """ Setting the parameter/coefficient to the Structure. The parameters to be set are given as keyword arguments.

                Args:
                    **kwargs: The parameters which are going to be set. Accepted parameters are: 'linker' ('dna'/None),
                    'source' ('experiment'/'cg'/'aa'), 'speed' (float), 'states' (int), 'residues_distance', (float),
                    'low_force_cutoff' (float), 'initial_guess' (dictionary), 'plot_columns' (int).

                Returns:
                    True if successful, False otherwise.
                """

        self.parameters = {**self.parameters, **kwargs}
        if overide:
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
        # self.make_partial_plots()
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
        max_contour_length = self.parameters['boundaries'][1]  # to align plots

        for n in range(1 + rows//5):
            r = min(5, rows-5*n)
            if r == 0:
                continue
            fig = plt.subplots(dpi=300, figsize=(10*int(columns), 5*r))
            axes = []

            # plotting each trace
            for k in range(5*n*columns, min(number, 5*columns*(n+1))):
                nk = k - 5*n*columns
                axes.append(plt.subplot2grid((r, 2 * columns), (int(nk / columns), (2 * nk) % (2 * columns))))
                self.traces[k].plot_histogram(position=axes[-1], max_contour_length=max_contour_length)
                axes.append(plt.subplot2grid((r, 2 * columns), (int(nk / columns), ((2 * nk) + 1) % (2 * columns))))
                self.traces[k].plot_fits(position=axes[-1])
            plt.tight_layout()

            if not output:
                ofile = self.name + '_contour_lengths_' + str(n+1) +'.png'
            else:
                ofile = output

            if self.debug:
                logger = set_logger(self.name)
                logger.info("Saving fits figure to " + ofile)
                close_logs(logger)

            plt.savefig(ofile)
            plt.close()
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

        fig = plt.figure(dpi=600, figsize=(10, 10))
        gs = fig.add_gridspec(4, 4)
        axes = [fig.add_subplot(gs[:2, :2]),
                fig.add_subplot(gs[:2, 2:]),
                fig.add_subplot(gs[2, :2]),
                fig.add_subplot(gs[3, :2]),
                fig.add_subplot(gs[2:, 2:])]

        #axes = [fig.add_subplot(gs[:2, :2]),
        #        fig.add_subplot(gs[:2, 2:]),
        #        fig.add_subplot(gs[2:, :2]),
        #        fig.add_subplot(gs[2:, 2:])]

        self._plot_total_contour_length_histo(axes[0])   # the total contour length histogram
        self._plot_overlaid_traces(axes[1])              # the overlaid smoothed traces
        self._plot_contour_length_histo(axes[2])
        self._plot_forces_histogram(axes[3])             # the rupture forces histogram
        #self._plot_dhs_analysis(axes[4])                 # Dudko analysis
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
        self._plot_dhs_analysis(axes)

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
        if os.path.isfile(self.parameters.get('state', '.') + '/' + self.name + '_hist_values'):
            with open(self.parameters.get('state', '.') + '/' + self.name + '_hist_values', 'r') as myfile:
                d = myfile.read()
            d = d.strip('[]\n')
            self.hist_values += [float(x) for x in d.split(', ')]
        else:
            for t in self.traces:
                self.hist_values += list(t.data['hist_values'])

            with open(self.parameters.get('state', '.') + '/' + self.name + '_hist_values', 'w') as myfile:
                myfile.write(str(self.hist_values))

            print("written hist_values")

        # decomposing the contour length histogram = obtaining the contour lengths
        max_length = (self.parameters['residues'] - 1) * self.parameters['residues_distance'] + 20
        print("collecting coefficients")
        parameters, boundaries = decompose_histogram(np.array(self.hist_values), states=self.num_states,
                                                     significance=self.parameters['significance'], max_value=max_length,
                                                     background_level=0.005, init_means=self.parameters['init_means'])
        file = open('cl_histo.txt', 'a')
        file.write(str(parameters) + '\n')
        file.close()

        self.parameters['l_prot'] = parameters
        self.parameters['boundaries'] = boundaries

        # collecting the rupture forces
        self._analyze_rupture_forces()

        # perform Dudko-Hummer-Szabo analysis
        # self._analyze_dhs()
        #
        # self.plot_dhs()

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
                'Residue distance (nm)\t' + str(self.parameters['residues_distance']),
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
    def _read_data(self, input_data=None, cases=None, columns=None, trace_name=None, **kwargs):
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
            if trace_name == None and headers[k][0] == 'd':
                trace_name = self.name + '_' + headers[k].strip('d_')
            elif trace_name == None:
                trace_name = self.name + '_' + str(int(k/2))
            else:
                pass
            trace_data = data[[headers[k], headers[k + 1]]]
            new_headers = {headers[k]: 'd', headers[k + 1]: 'F'}
            trace_data = trace_data.rename(columns=new_headers)
            if parameters['unit'] == 'A':
                trace_data['d'] = trace_data['d'] / 10
            trace_data = trace_data[trace_data['F'] < self.parameters['max_force']]
            self.traces.append(Trace(trace_name, trace_data, input_data, experiment_name=self.name, debug=self.debug,
                                     parameters=parameters))
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
        self.parameters = deepcopy(default_parameters)

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
        self.num_states = None
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
        data_to_plot = [x for x in self.hist_values if x > 0]
        bins = 4 * int(max(data_to_plot))
        position.hist(data_to_plot, bins=bins, density=True, alpha=0.5)

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
        max_f, min_d, max_d = 0, np.inf, 0

        # plotting overlaid smoothed traces
        for k in range(len(self.traces)):
            t = self.traces[k]
            position.plot(t.smoothed['d'], t.smoothed['F'], color=get_gray_shade(k, len(self.traces)))
            max_f = max(max_f, t.data['F'].max())
            min_d = min(min_d, t.data['d'].min())
            max_d = max(max_d, t.data['d'].max())

        # plotting fits
        plot_trace_fits(position, self.parameters, max_f, self.parameters['residues_distance'],
                        method=self.parameters['method'])

        position.set_ylim(0, max_f)
        position.set_xlim(min_d, max_d)
        position.legend()
        return

    def _plot_contour_length_histo(self, position):
        if len(self.forces) == 0:
            if self.debug:
                logger = set_logger(self.name)
                logger.debug("Structure.forces empty, nothing to plot. Did you run .analyze() or "
                             ".collect_coefficients ?")
                close_logs(logger)
            return

        if self.debug:
            logger = set_logger(self.name)
            logger.info("Making the contour length histogram.")
            close_logs(logger)

        if self.parameters['boundaries']:
            bound = self.parameters['boundaries'][1]
        else:
            bound = max(self.hist_values)

        position.set_xlabel('Contour length [nm]')
        position.set_ylabel('Occurences')
        position.set_title('Histogram of contour lengths')
        position.set_xlim(0, bound)
        # position.set_ylim(0, 2)

        lspace = np.linspace(0, bound)

        k = 0
        print("plotting contour length histo")
        all_data = []
        for state in self.states:
            print(state)
            data_to_plot = list(self.forces[self.forces['state'] == state]['means'].dropna())
            all_data = data_to_plot + all_data

            if len(data_to_plot) < 3:
                continue
        all_data = np.array(all_data)
        parameters, boundaries = decompose_histogram(all_data, significance=self.parameters['significance'],
                                                     states=1)

        print(parameters)
        # if len(parameters) == 0:
        #    continue

        mean, width, height = parameters[['means', 'widths', 'heights']].values[0]
        #label = 'Mean: ' + str(round(mean, 3)) + ' nm'

        bins = max(int(max(all_data)) - int(min(all_data)), 1)
        position.hist(all_data, bins=bins, alpha=0.5, density=True)

        for index, row in parameters[['means', 'widths', 'heights']].iterrows():
            mean, width, height = tuple(row.to_numpy())
            label = "Mean: " + str(round(mean, 3))
            y_plot = single_gaussian(lspace, height, mean, width)
            position.plot(lspace, y_plot, linestyle='--', linewidth=0.5, label=label,
                          color=get_color(k, len(parameters)))
            k += 1




        position.legend()
        return

    def _plot_forces_histogram(self, position):
        if len(self.forces) == 0:
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

        print(self.parameters['l_prot'])

        # plotting the histograms and the fitted contour
        k = 0
        print("plotting force histogram")
        for state in self.states:
            print(state)
            data_to_plot = np.array(self.forces[self.forces['state'] == state]['rupture_forces'].dropna())
            if len(data_to_plot) < 3:
                continue

            file = open('rupture.txt', 'a')
            file.write('State ' + str(state) + '\n')
            file.write(str(data_to_plot) + '\n')
            file.close()

            parameters, boundaries = decompose_histogram(data_to_plot, significance=self.parameters['significance'],
                                                         states=1)
            print(parameters)
            file = open('rupture_gauss.txt', 'a')
            file.write(str(state) + '\n')
            file.write(str(parameters) + '\n')
            file.close()

            mean, width, height = parameters[['means', 'widths', 'heights']].values[0]
            label = str(round(mean, 3)) + ' pN'

            bins = max(int(max(data_to_plot)) - int(min(data_to_plot)), 1)
            position.hist(data_to_plot, bins=bins, color=get_color(k, len(self.states)), alpha=0.5, density=True,
                          label=label)
            y_plot = single_gaussian(f_space, height, mean, width)
            position.plot(f_space, y_plot, linestyle='--', linewidth=0.5, color=get_color(k, len(self.states)))
            k += 1
        position.legend()
        return

    def _plot_dhs_analysis(self, position):
        # setting the scene
        position.set_title('Dudko-Hummer-Szabo lifetime')
        position.set_xlabel('Rupture force [pN]')
        position.set_ylabel('State lifetime [s]')
        # position.set_ylim(top=1000)
        position.set_yscale('log')

        for k, data in enumerate(self.dhs_states):
            beg, end, l_prot, dhs_data = data
            if len(dhs_data) == 0:
                continue
            label = 'l_prot=' + str(round(l_prot, 3))
            position.plot(dhs_data['forces'], dhs_data['lifetime'], label=label,
                          color=get_color(k, len(self.dhs_states)))
        position.legend()
        return

    # analyzing
    def _analyze_rupture_forces(self):
        """ Collecting and analyzing rupture forces"""

        # collecting forces
        print("collecting forces")
        self.forces = pd.concat([t.parameters['l_prot'][['means', 'rupture_forces']].dropna() for t in self.traces],
                                    ignore_index=True)
        self.max_f = self.forces['rupture_forces'].max()

        # assigning a rupture force to the state from the list of cumulative states
        for index, row in self.parameters['l_prot'][['means', 'widths', 'heights']].iterrows():
            mean, width, height = tuple(row.to_numpy())
            self.forces['state_' + str(index)] = single_gaussian(np.array(self.forces['means']), height, mean, width)
            self.states.append('state_' + str(index))

        # selecting the corresponding state
        self.forces['state'] = self.forces[self.states].idxmax(axis=1)
        print(self.forces)
        return True

    def _analyze_dhs(self):
        # TODO clean it up
        """ Calculating the data for the Dudko-Hummer-Szabo analysis of the states lifetime"""
        print("analyzing dhs")
        for k in range(len(self.states)):
            state = self.states[k]
            print(state)

            # obtaining the force distribution
            data_to_plot = np.array(self.forces[self.forces['state'] == state]['rupture_forces'].dropna())
            if len(data_to_plot) < 2:
                continue
            parameters, boundaries = decompose_histogram(data_to_plot, significance=self.parameters['significance'],
                                                         states=1)
            print(parameters)
            mean, width, height = parameters[['means', 'widths', 'heights']].values[0]
            f_space = np.linspace(*boundaries)

            # finding protein contour length
            state_index = int(state.strip('state_'))
            l_prot = self.parameters['l_prot'].loc[state_index, 'means']

            # data for dhs
            dhs_data = pd.DataFrame({'forces': f_space,
                                     'force_load': get_force_load(f_space, self.parameters),
                                     'probability': single_gaussian(f_space, height, mean, width),
                                     'nominator': integrate_gauss(f_space, mean, width)})
            dhs_data['denominator'] = dhs_data['probability'] * dhs_data['force_load']
            dhs_data = dhs_data[dhs_data['denominator'] > 0.1]
            dhs_data['lifetime'] = dhs_data['nominator'] / dhs_data['denominator']
            if len(dhs_data) < 4:
                self.dhs_states.append(boundaries + [l_prot, []])
                continue

            print(dhs_data)

            # calculating the lifetime and fitting it
            self.dhs_states.append(boundaries + [l_prot, dhs_data])
            self.dhs_results[round(l_prot, 3)] = dhs_feat(dhs_data, l_prot)
            print(self.dhs_results[round(l_prot, 3)])
        return True
