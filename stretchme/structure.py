from .stretching_tools import *
from .default_parameters import default_parameters
from .trace import Trace
from matplotlib import pyplot as plt
import logging
import os
from io import StringIO


class Structure:
    """The class containing information about the whole experiment, consisting of many measurements (traces).

    Attributes:
        input_data (str/Pandas Dataframe/None): The data to be analyzed or the path the file with data.
        parameters (dict): The dictionary of measurement parameters.
        coefficients (dict): The dictionary with fitted coefficients.
        rupture_forces (dict): The dictionary with measured rupture forces.

    """
    def __init__(self, input_data=None, cases=None, columns=None, parameters={}, name=None, debug=False, **kwargs):
        self.orig_input = input_data
        # setting the name
        self._set_name(input_data, name)

        # initializing log file if debug
        self.logger = None
        if debug:
            logging.basicConfig(filename=self.name + '.log', filemode='w', level=logging.DEBUG,
                                format='%(asctime)s %(message)s')
            self.logger = logging.getLogger()
            self.logger.debug('Initializing the "Structure" class with name: ' + str(self.name))

        # setting parameters
        self._set_parameters(parameters, **kwargs)

        # reading the data (setting the traces)
        self.traces = []
        self._read_data(input_data, cases, columns)
        return

    def _close_logs(self):
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)
        return

    def _read_data(self, input_data, cases=None, columns=None, **kwargs):
        # TODO clean up
        if 'separator' in kwargs.keys():
            separator = kwargs['separator']
            del kwargs['separator']
        else:
            separator = self.parameters['separator']
        if 'sheet_name' in kwargs.keys():
            sheet_name = kwargs['sheet_name']
            del kwargs['sheet_name']
        else:
            sheet_name = self.parameters['sheet_name']
        if (isinstance(input_data, pd.DataFrame) and len(input_data) == 0) or \
                (not isinstance(input_data, pd.DataFrame) and not input_data):
            if self.logger:
                self.logger.debug("Initializing empty class. Hope you'll add some traces to analyze.")
            return

        # if there is the input_data to read
        if isinstance(input_data, pd.DataFrame):
            data = input_data
            filename = 'Pandas Dataframe'
        elif isinstance(input_data, str) and os.path.isfile(input_data):
            data = self._read_input_file(input_data, separator, sheet_name, cases, columns, **kwargs)
            filename = input_data
        else:
            raise NotImplementedError("Could not understand the input. Currently reading Pandas Dataframe, and paths "
                                      "to .xls and .csv files.")
        data.dropna(axis='columns', how='all', inplace=True)
        headers = list(data)
        if columns:
            if isinstance(columns[0], int) and isinstance(columns[0], int) \
                    and (columns[0] not in headers or columns[1] not in headers):
                trace_data = data[[headers[columns[0]], headers[columns[1]]]]
                new_headers = {headers[columns[0]]: 'd', headers[columns[1]]: 'F'}
            elif columns[0] in headers and columns[1] in headers:
                trace_data = data[columns]
                new_headers = {columns[0]: 'd', columns[1]: 'F'}
            else:
                raise NotImplementedError("Could not understand the columns to be read.")

            trace_name = len(self.traces)
            trace_data = trace_data.rename(columns=new_headers)
            if 'unit' in kwargs.keys() and kwargs['unit'] == 'A':
                trace_data['d'] = trace_data['d'] / 10
            self.traces.append(Trace(trace_name, trace_data, filename, logger=self.logger,
                                     parameters=self._set_trace_parameters(kwargs)))
        else:
            for k in range(len(headers)%2, len(headers), 2):
                # TODO take care of "index list out of range here, if the data are given in wrong format"
                trace_data = data[[headers[k], headers[k + 1]]]
                new_headers = {headers[k]: 'd', headers[k + 1]: 'F'}
                trace_data = trace_data.rename(columns=new_headers)
                if headers[k][0] == 'd':
                    trace_name = headers[k].strip('dF_ ')
                else:
                    trace_name = len(self.traces)
                if 'unit' in kwargs.keys() and kwargs['unit'] == 'A':
                    trace_data['d'] = trace_data['d'] / 10
                if not cases or trace_name in cases:
                    self.traces.append(Trace(trace_name, trace_data, filename, logger=self.logger,
                                     parameters=self._set_trace_parameters(kwargs)))
        return

    def _read_input_file(self, filename, separator, sheet_name, cases, columns, **kwargs):
        if self.logger:
            self.logger.info(
                'Reading ' + filename + ' cases ' + str(cases) + ' columns ' + str(columns) + str(kwargs))
        if '.xls' in filename:
            data = pd.read_excel(filename, sheet_name=sheet_name)
        else:  # the .csv-like case
            if '.csv' not in filename:
                if self.logger:
                    self.logger.warning(
                        "Treating the input file as .csv file, although it does not have such extension.")
            with open(filename, 'r') as myfile:
                content = myfile.read().split("#")[1].strip()
            if separator != ' ':
                data = pd.read_csv(StringIO(content), sep=separator, escapechar='#')
            else:
                data = pd.read_csv(StringIO(content), delim_whitespace=True, escapechar='#')
        return data

    def add_trace(self, filename, cases=None, columns=None, **kwargs):
        self._read_data(filename, cases, columns, **kwargs)
        return

    # setting parameters
    def _set_name(self, filename, name):
        if name:
            self.name = name
        elif isinstance(filename, str) and len(filename) > 0:
            self.name = os.path.splitext(os.path.basename(filename))[0]
        else:
            self.name = 'Experiment'
        return 0

    def _set_parameters(self, parameters, **kwargs):
        self.parameters = default_parameters
        for key in parameters.keys():
            self.parameters[key] = parameters[key]
        for key in kwargs:
            self.parameters[key] = kwargs[key]
        if self.logger:
            self.logger.debug("Parameters:\t" + str(parameters))
        return

    def _set_trace_parameters(self, additional={}):
        # TODO clean up
        parameters = {}
        for key in ['linker', 'source', 'speed', 'residues_distance', 'states', 'low_force_cutoff']:
            parameters[key] = self.parameters[key]
            if key in additional.keys():
                parameters[key] = additional[key]
        key = 'initial_guess'
        if key in additional.keys():
            parameters[key] = additional[key]
        elif key in self.parameters.keys() and isinstance(self.parameters[key], dict):
            parameters[key] = self.parameters[key]
        else:
            parameters[key] = default_parameters[key][parameters['source']]
        return parameters

    def set_residues(self, value):
        self.parameters['residues'] = value
        if self.logger:
            self.logger.info("Set 'residues' value to " + str(value))
        return

    def set_distance(self, value):
        self.parameters['distance'] = value
        if self.logger:
            self.logger.info("Set 'distance' value to " + str(value))
        return

    def set_linker(self, value):
        self.parameters['linker'] = value
        if self.logger:
            self.logger.info("Set 'linker' value to " + str(value))
        return

    def set_source(self, value):
        self.parameters['source'] = value
        if self.logger:
            self.logger.info("Set 'source' value to " + str(value))
        return

    def set_unit(self, value):
        self.parameters['unit'] = value
        if self.logger:
            self.logger.info("Set 'unit' value to " + str(value))
        return

    def set_speed(self, value):
        self.parameters['speed'] = value
        if self.logger:
            self.logger.info("Set 'speed' value to " + str(value))
        return

    def set_residues_distance(self, value):
        self.parameters['residues_distance'] = value
        if self.logger:
            self.logger.info("Set 'residues_distance' value to " + str(value))
        return

    def set_minimal_stretch_distance(self, value):
        self.parameters['minimal_stretch_distance'] = value
        if self.logger:
            self.logger.info("Set 'minimal_stretch_distance' value to " + str(value))
        return

    def set_plot_columns(self, value):
        self.parameters['plot_columns'] = value
        if self.logger:
            self.logger.info("Set 'plot_columns' value to " + str(value))
        return

    def set_high_force_cutoff(self, value):
        self.parameters['high_force_cutoff'] = value
        if self.logger:
            self.logger.info("Set 'high_force_cutoff' value to " + str(value))
        return

    def set_low_force_cutoff(self, value):
        self.parameters['low_force_cutoff'] = value
        if self.logger:
            self.logger.info("Set 'low_force_cutoff' value to " + str(value))
        return

    def set_max_rupture_force(self, value):
        self.parameters['max_rupture_force'] = value
        if self.logger:
            self.logger.info("Set 'max_rupture_force' value to " + str(value))
        return

    def set_initial_guess(self, value):
        self.parameters['initial_guess'] = value
        if self.logger:
            self.logger.info("Set 'initial_guess' value to " + str(value))
        return

    def set_states(self, value):
        self.parameters['states'] = value
        if self.logger:
            self.logger.info("Set 'states' value to " + str(value))
        return

    # analyzing
    def _collect_coefficients(self):
        self.coefficients = {}
        self.coefficients['p_prot'] = np.mean(np.array([t.coefficients.get('p_prot', 0) for t in self.traces]))
        self.coefficients['p_dna'] = np.mean(np.array([t.coefficients.get('p_dna', 0) for t in self.traces]))
        self.coefficients['l_dna'] = np.mean(np.array([t.coefficients.get('l_dna', 0) for t in self.traces]))
        k_prots = [t.coefficients.get('k_prot', None) for t in self.traces]
        if None in k_prots:
            self.coefficients['k_prot'] = None
        else:
            self.coefficients['k_prot'] = np.mean(np.array(k_prots))
        k_dnas = [t.coefficients.get('k_dna', None) for t in self.traces]
        if None in k_dnas:
            self.coefficients['k_dna'] = None
        else:
            self.coefficients['k_dna'] = np.mean(np.array(k_dnas))
        return

    def _find_paths(self):
        return

    def analyze(self):
        if self.logger:
            self.logger.info("Starting analysis of the experiment.")
        for t in self.traces:
            t.analyze()
        self._collect_coefficients()
        self.make_plots()
        self._find_paths()
        # self.save_data()
        return

    # plotting
    def make_plots(self):
        self._make_histograms()
        self._make_partial_plots()
        return

    def _make_partial_plots(self):
        if self.logger:
            print("Making plots of individual fits.")
        # setting the figure
        number = len(self.traces)                                   # number of traces to plot
        columns = self.parameters['plot_columns']                   # number of columns
        rows = max(int(np.ceil(float(number) / columns)), 2)        # number of rows
        fig = plt.subplots(dpi=600, figsize=(10*int(columns), 5*int(rows)))
        axes = []
        max_contour_length = self.coefficients['boundaries'][1]     # to align plots

        # plotting each trace
        for k in range(0, number):
            axes.append(plt.subplot2grid((rows, 2 * columns), (int(k / columns), (2 * k) % (2 * columns))))
            self.traces[k]._plot_histogram(position=axes[-1], max_contour_length=max_contour_length)
            axes.append(plt.subplot2grid((rows, 2 * columns), (int(k / columns), ((2 * k) + 1) % (2 * columns))))
            self.traces[k]._plot_fits(position=axes[-1])

        plt.tight_layout()

        if self.logger:
            print("Saving contour lengths figure to " + self.name + '_contour_lengths.png')
        plt.savefig(self.name + '_contour_lengths.png')
        return

    def _plot_total_contour_length_histo(self, position, significance=0.001):
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')

        # collecting the data from all traces
        hist_values = []
        for t in self.traces:
            hist_values += list(t.data['hist_values'])

        # plotting histogram
        position.hist(hist_values, bins=500, density=True, alpha=0.5)

        # decomposing histogram
        parameters, boundaries = decompose_histogram(np.array(hist_values), significance, self.parameters['states'])
        self.coefficients['l_prot'] = parameters
        self.coefficients['boundaries'] = boundaries
        position.set_xlim(0, boundaries[1])

        # plotting decomposed histogram
        l_space = np.linspace(0, boundaries[1], 1001)
        plot_decomposed_histogram(position, self.coefficients['l_prot'], l_space, self.parameters['residues_distance'])

        position.legend()
        return

    def _plot_overlaid_traces(self, position):
        position.set_xlabel('Distension [nm]')
        position.set_ylabel('Force [pN]')
        position.set_title('All smoothed traces overlaid')
        max_f = 0

        # plotting overlaid smoothed traces
        for k in range(len(self.traces)):
            t = self.traces[k]
            position.plot(t.smoothed['d'], t.smoothed['F'], color=get_gray_shade(k, len(self.traces)))
            max_f = max(max_f, t.data['F'].max())

        position.set_ylim(0, max_f)

        # plotting fits
        f_space = np.linspace(0.1, max_f)
        plot_trace_fits(position, self.coefficients, f_space, self.parameters['residues_distance'])

        position.legend()
        return

    def _plot_forces_histogram(self, position, significance=0.001):
        position.set_xlabel('Force [pN]')
        position.set_ylabel('Occurences')
        position.set_title('Rupture force histogram')

        # collecting forces with respective contour lengths
        forces = pd.concat([t.coefficients['l_prot'][['means', 'rupture_forces']].dropna() for t in self.traces],
                           ignore_index=True)
        self.max_f = forces['rupture_forces'].max()
        columns = ['l_prot', 'p_prot', 'k_prot', 'means', 'widths', 'heights', 'beg', 'end']
        self.forces_cofficients = pd.DataFrame(columns=columns)
        f_space = np.linspace(0, self.max_f)
        states = []

        # calculating the probability of force contour length coming from calculated contour length for whole experiment
        for index, row in self.coefficients['l_prot'][['means', 'widths', 'heights']].iterrows():
            mean, width, height = tuple(row.to_numpy())
            forces['state_' + str(index)] = gauss(np.array(forces['means']), height, mean, width)
            states.append('state_' + str(index))

        # selecting the corresponding state
        forces['state'] = forces[states].idxmax(axis=1)

        # plotting the histogram and the fitted contour
        for k in range(len(states)):
            state = states[k]
            data = np.array(forces[forces['state'] == state]['rupture_forces'].dropna())
            parameters, boundaries = decompose_histogram(data, significance)
            y_plot = np.zeros(len(f_space))
            for index, row in parameters.iterrows():
                y_plot += gauss(f_space, row['heights'], row['means'], row['widths'])
            label = '; '.join([str(round(x, 3)) for x in list(parameters['means'].values)]) + ' pN'
            position.hist(data, bins=20, color=get_color(k, len(states)), alpha=0.5, density=True, label=label)
            position.plot(f_space, y_plot, linestyle='--', color=get_color(k, len(states)))

            # preparing data for Dudko analysis
            state_index = int(state.strip('state_'))
            l_prot = self.coefficients['l_prot'].loc[state_index, 'means']
            parameters['l_prot'] = np.array([l_prot for _ in range(len(parameters))])
            parameters['p_prot'] = np.array([self.coefficients['p_prot'] for _ in range(len(parameters))])
            parameters['k_prot'] = np.array([self.coefficients['k_prot'] for _ in range(len(parameters))])
            parameters['beg'] = np.array([boundaries[0] for _ in range(len(parameters))])
            parameters['end'] = np.array([boundaries[1] for _ in range(len(parameters))])
            self.forces_cofficients = pd.concat([self.forces_cofficients, parameters], ignore_index=True)

        position.legend()
        return

    def _plot_dhs_analysis(self, position):
        self.dhs_results = {}
        position.set_title('Dudko-Hummer-Szabo lifetime')
        position.set_xlabel('Rupture force [pN]')
        position.set_ylabel('State lifetime [s]')
        position.set_yscale('log')
        if not self.parameters['speed']:
            self.parameters['speed'] = 1
            if self.logger:
                self.logger.info("Unknown extension speed. Set extension speed to 1. "
                                 "You may set the speed using Experiment.set_speed() function")

        k = 0
        columns = ['l_prot', 'p_prot', 'k_prot', 'means', 'widths', 'heights', 'beg', 'end']
        for index, row in self.forces_cofficients[columns].iterrows():
            l_prot, p_prot, k_prot, mean, width, height, beg, end = tuple(row.to_numpy())
            f_space = np.linspace(beg, end, 100)

            force_load = get_force_load(f_space, self.parameters['speed'], l_prot, p_prot, k_prot)
            probability = gauss(f_space, height, mean, width)
            denominator = probability * force_load
            nominator = integrate_gauss(f_space, mean, width)

            # calculating the lifetime and fitting it
            lifetime = nominator / denominator

            f_space = f_space[lifetime < 10000]
            lifetime = lifetime[lifetime < 10000]
            parameters = dhs_feat(f_space, lifetime)
            self.dhs_results[(round(float(row['l_prot']), 3), round(float(row['means']), 3))] = parameters
            label = 'l_prot=' + str(round(float(row['l_prot']), 3)) + '; force=' + str(round(float(row['means']), 3))
            position.plot(f_space, lifetime, label=label, color=get_color(k, len(self.forces_cofficients)))
            k += 1

        position.legend()
        return

    def _make_histograms(self):
        if self.logger:
            self.logger.info("Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        self._plot_total_contour_length_histo(axes[0, 0])   # the total contour length histogram
        self._plot_overlaid_traces(axes[0, 1])              # the overlaid smoothed traces
        # TODO take care of extremal cases with only one fit for the forces and Dudko
        try:
            self._plot_forces_histogram(axes[1, 0])             # the rupture forces histogram
        except:
            pass
        try:
            self._plot_dhs_analysis(axes[1, 1])               # Dudko analysis
        except:
            pass

        fig.tight_layout()
        if self.logger:
            self.logger.info("Saving histograms figure to " + self.name + '_histograms.png')
        plt.savefig(self.name + '_histograms.png')
        return

    # printing the data
    def save_data(self):
        oname = str(self.name) + '_results'
        separator = '################\n'

        if self.logger:
            self.logger.info("Saving output to " + oname)
        result = []

        # general info
        result.append('General info:')
        result.append('Name:\t\t\t' + str(self.name))
        result.append('Residues:\t\t\t' + str(self.parameters['residues']))
        result.append('Linker:\t\t\t' + str(self.parameters['linker']))
        result.append('Unit:\t\t\t' + str(self.parameters['unit']))
        result.append('Data source:\t\t' + str(self.parameters['source']))
        result.append('Pulling speed:\t\t\t' + str(self.parameters['speed']))
        result.append(separator)

        # summary of individual curve
        result.append('Summary of individual curves')
        for t in self.traces:
            result.append(t.summary())
        result.append(separator)

        # summary of the cumulative statistics
        result.append('Summary of the cumulative statistics')
        result.append('p_Prot:\t\t' + str(self.coefficients['p_prot']))
        result.append('p_DNA:\t\t' + str(self.coefficients['p_dna']))
        result.append('k_Prot:\t\t' + str(self.coefficients['k_prot']))
        result.append('k_DNA:\t\t' + str(self.coefficients['k_dna']))
        result.append('l_DNA:\t\t' + str(self.coefficients['l_dna']))
        result.append('Contour length\tgamma\tsigma^2')
        result.append(self.coefficients['l_prot'].to_csv(sep='\t'))
        result.append('Contour length boundaries')
        result.append(str(self.coefficients['boundaries']))
        result.append(separator)

        # Dudko-Hummer-Szabo analysis
        result.append('Dudko-Hummer-Szabo analysis')
        result.append(self.forces_cofficients.to_csv(sep='\t'))
        for key in self.dhs_results.keys():
            for v in self.dhs_results[key].keys():
                result.append(str(key) + '\t' + str(v) + '\t' + str(self.dhs_results[key][v]))
        result.append(separator)

        with open(oname, 'w') as ofile:
            ofile.write('\n'.join(result))
        # TODO make something intelligent with the loggers, to close it without the need to run 'save_data'
        self._close_logs()
        return

    def get_info(self):
        info = []
        info.append('Experiment name ' + str(self.name))
        info.append('Experiment source file ' + str(self.orig_input))
        info.append('Number of traces: ' + str(len(self.traces)))
        info.append('Details of traces: ' + str(len(self.traces)))
        for t in self.traces:
            info.append(t.get_info())
        return '\n'.join(info)
