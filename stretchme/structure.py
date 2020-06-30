from .stretching_tools import *
from .default_parameters import default_parameters
from .trace import Trace
from matplotlib import pyplot as plt
import logging
import os
from io import StringIO


class Structure:
    def __init__(self, filename=None, cases=None, columns=None, parameters={}, name=None, debug=False):
        self.orig_input = filename
        # setting the name
        self._set_name(filename, name)

        # initializing log file if debug
        self.logger = None
        if debug:
            logging.basicConfig(filename=self.name + '.log', filemode='w', level=logging.DEBUG,
                                format='%(asctime)s %(message)s')
            self.logger = logging.getLogger()
            self.logger.debug('Initializing the "Structure" class with name: ' + str(self.name))

        # setting parameters
        self._set_parameters(parameters)

        # reading the data (setting the traces)
        self.traces = []
        self._read_data(filename, cases, columns)
        return

    def _read_data(self, filename, cases=None, columns=None, **kwargs):
        # TODO clean up
        if 'separator' in kwargs.keys():
            separator = kwargs['separator']
        else:
            separator = self.parameters['separator']
        if 'sheet_name' in kwargs.keys():
            sheet_name = kwargs['sheet_name']
        else:
            sheet_name = self.parameters['sheet_name']
        if not filename:
            if self.logger:
                self.logger.debug("Initializing empty class. Hope you'll add some traces to analyze.")
            return

        # if there is a file to read
        if self.logger:
            self.logger.info('Reading ' + filename + ' cases ' + str(cases) + ' columns ' + str(columns) + str(kwargs))
        if '.xls' in filename:
            data = pd.read_excel(filename, sheet_name=sheet_name)
        else:   # the .csv-like case
            if '.csv' not in filename:
                if self.logger:
                    self.logger.warning("Treating the input file as .csv file, although it does not have such extension.")
            with open(filename, 'r') as myfile:
                content = myfile.read().split("#")[1].strip()
            if separator != ' ':
                data = pd.read_csv(StringIO(content), sep=separator, escapechar='#')
            else:
                data = pd.read_csv(StringIO(content), delim_whitespace=True, escapechar='#')
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
            for k in range(0, len(headers), 2):
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
                if cases and trace_name in cases:
                    self.traces.append(Trace(trace_name, trace_data, filename, logger=self.logger,
                                     parameters=self._prepare_trace_parameters(kwargs)))
        return

    def add_trace(self, filename, cases=None, columns=None, **kwargs):
        self._read_data(filename, cases, columns, **kwargs)
        return

    # setting parameters
    def _set_name(self, filename, name):
        if name:
            self.name = name
        elif filename:
            self.name = os.path.splitext(os.path.basename(filename))[0]
        else:
            self.name = 'Experiment'
        return 0

    def _set_parameters(self, parameters):
        self.parameters = default_parameters
        for key in parameters.keys():
            self.parameters[key] = parameters[key]
        if self.logger:
            self.logger.debug("Parameters:\t" + str(parameters))
        return

    def _set_trace_parameters(self, additional={}):
        # TODO clean up
        parameters = {}
        for key in ['residues', 'distance', 'linker', 'source', 'speed', 'residues_distance',
                    'minimal_stretch_distance', 'max_cluster_gap', 'initial_guess']:
            parameters[key] = self.parameters[key]
            if key in additional.keys():
                parameters[key] = additional[key]
        for key in ['high_force_cutoff', 'low_force_cutoff', 'max_rupture_force']:
            if key in additional.keys():
                parameters[key] = additional[key]
            elif key in self.parameters.keys() and isinstance(self.parameters[key], str):
                parameters[key] = self.parameters[key]
            elif key in self.parameters.keys() and isinstance(self.parameters[key], dict):
                parameters[key] = self.parameters[key][parameters['source']]
            else:
                parameters[key] = default_parameters[key][parameters['source']]
        return parameters

    def set_residues(self, value):
        self.parameters['residues'] = value
        return

    def set_distance(self, value):
        self.parameters['distance'] = value
        return

    def set_linker(self, value):
        self.parameters['linker'] = value
        return

    def set_source(self, value):
        self.parameters['source'] = value
        return

    def set_unit(self, value):
        self.parameters['unit'] = value
        return

    def set_speed(self, value):
        self.parameters['speed'] = value
        return

    def set_residues_distance(self, value):
        self.parameters['residues_distance'] = value
        return

    def set_minimal_stretch_distance(self, value):
        self.parameters['minimal_stretch_distance'] = value
        return

    def set_max_cluster_gap(self, value):
        self.parameters['max_cluster_gap'] = value
        return

    def set_plot_columns(self, value):
        self.parameters['plot_columns'] = value
        return

    def set_high_force_cutoff(self, value):
        self.parameters['high_force_cutoff'] = value
        return

    def set_low_force_cutoff(self, value):
        self.parameters['low_force_cutoff'] = value
        return

    def set_max_rupture_force(self, value):
        self.parameters['max_rupture_force'] = value
        return

    def set_initial_guess(self, value):
        self.parameters['initial_guess'] = value
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
        number = len(self.traces)       # number of traces to plot
        columns = self.parameters['plot_columns']
        rows = max(int(np.ceil(float(number) / columns)), 2)
        fig, axes = plt.subplots(rows, 2 * columns, dpi=600, figsize=(10*int(columns), 5*int(rows)))
        max_contour_length = self.coefficients['boundaries'][1]
        for k in range(0, number):
            self.traces[k]._plot_histogram(position=axes[int(k / columns), (2 * k) % (2 * columns)],
                                           max_contour_length=max_contour_length)
            self.traces[k]._plot_fits(position=axes[int(k / columns), ((2 * k) + 1) % (2 * columns)])
        fig.tight_layout()
        if self.logger:
            print("Saving contour lengths figure to " + self.name + '_contour_lengths.png')
        plt.savefig(self.name + '_contour_lengths.png')
        return

    def _plot_individual_contour_length_histo(self, position):
        position.set_xlabel('Fitted contour length')
        position.set_ylabel('Occurences')
        position.set_title('Histogram of fitted contour lengths')
        lengths = [parameter[0] for t in self.traces for parameter in t.coefficients['lengths']]
        position.hist(lengths, bins=50)
        return

    def _plot_total_contour_length_histo(self, position, boundary_value=0.001):
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')

        # collecting the data from all traces
        hist_values = []
        for t in self.traces:
            hist_values += list(t.data['hist_values'])

        # plotting histogram
        position.hist(hist_values, bins=500, density=True)

        # decomposing histogram
        parameters, boundaries = decompose_histogram(np.array(hist_values), boundary_value)
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
        return

    def _plot_forces_histogram(self, position):
        position.set_xlabel('Force [pN]')
        position.set_ylabel('Occurences')
        position.set_title('Rupture force histogram')
        return

    def _plot_dudko_analysis(self, position):
        position.set_title('Dudko analysis')
        position.set_xlabel('Rupture force [pN]')
        position.set_ylabel('Rupture time [s]')
        # position.set_xlim(self.minf, self.maxf)
        # for mu in self.lengths:
        #     self.dudko_parameters[mu] = []
        #     heights = np.array(self.lengths[mu]['hist_vals'])
        #     df = self.lengths[mu]['hist_width']
        #     hist_beg = self.lengths[mu]['hist_beg']
        #     heights_div = np.array([heights[k] for k in range(len(heights)) if heights[k] != 0])
        #     x = np.array([hist_beg + (k - 1 / 2) * df for k in range(len(heights)) if heights[k] != 0])
        #     x_smooth = np.linspace(min(x), max(x))
        #     load = force_load(x, self.pprot, mu, self.parameters['speed'])
        #     area = np.array([heights[k - 1] / 2 + sum([heights[i - 1] for i in range(k + 1, len(heights)+1)])
        #                      for k in range(len(heights)) if heights[k] != 0])
        #     y = (area * df) / (load * heights_div)
        #     label = 'L=' + str(round(mu, 3))
        #     for v in self.parameters['Dudko parameters']:
        #         self.dudko_parameters[mu][v] = fit_dudko(x, y, v)
        #         t0 = self.dudko_parameters[mu][v]['t0']
        #         gamma = self.dudko_parameters[mu][v]['x'] / (2 * self.dudko_parameters[mu][v]['g'])
        #         g = self.dudko_parameters[mu][v]['g']
        #         y_smooth = t0 * (1 - gamma * x_smooth) ** (-1) * np.exp(-g * (1 - (1 - gamma * x_smooth) ** 2))
        #         position.plot(x, y, 'o', color=self.lengths[mu]['color'], label=label)
        #         position.plot(x_smooth, y_smooth, '-', color=self.lengths[mu]['color'])
        # position.legend()
        return

    def _make_histograms(self):
        if self.logger:
            print("Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        # the total contour length histogram
        self._plot_total_contour_length_histo(axes[0, 0])

        # the overlaid smoothed traces
        self._plot_overlaid_traces(axes[0, 1])

        # the rupture forces histogram
        self._plot_forces_histogram(axes[1, 0])

        # Dudko analysis
        self._plot_dudko_analysis(axes[1, 1])

        fig.tight_layout()
        if self.logger:
            print("Saving histograms figure to " + self.name + '_histograms.png')
        plt.savefig(self.name + '_histograms.png')
        return

    # printing the data
    def save_data(self):
        if self.logger:
            print("Saving output to " + self.name + '.log')
        result = []
        for t in self.traces:
            result.append(t.name)
            result.append(str(t.coefficients))
        with open(self.name + '.log', 'w') as ofile:
            ofile.write('\n'.join(result))
        # result = []
        # separator = '################\n'
        #
        # # general info
        # result.append('General info:')
        # result.append('Name:\t\t\t' + self.info['name'])
        # result.append('Residues:\t\t\t' + str(self.info['residues']))
        # result.append('End-to-end distance:\t' + str(self.info['distance']))
        # result.append('Linker:\t\t\t' + str(self.info['linker']))
        # result.append('Unit:\t\t\t' + str(self.info['unit']))
        # result.append('Data source:\t\t' + str(self.info['source']))
        # result.append('Pulling speed:\t\t\t' + str(self.info['speed']))
        # result.append(separator)
        #
        # # parameters
        # result.append('Calculation parameters:')
        # result.append('Data path:\t\t' + str(self.parameters['data_path']))
        # result.append('Data file prefix:\t\t' + str(self.parameters['data_file_prefix']))
        # result.append('Data file suffix:\t\t' + str(self.parameters['data_file_suffix']))
        # result.append('Residue-residue distance:\t' + str(self.parameters['residues_distance']))
        # result.append('Minimal distance between jumps:\t\t' + str(self.parameters['minimal_stretch_distance']))
        # result.append('Low force cutoff:\t\t' + str(self.parameters['low_force_cutoff']))
        # result.append('High force cutoff:\t\t' + str(self.parameters['high_force_cutoff']))
        # result.append('Minimal gap between peaks in cluster:\t\t' + str(self.parameters['cluster_max_gap']))
        # result.append('Number of columns in individual plots:\t\t' + str(self.parameters['columns']))
        # result.append(separator)
        #
        # # summary of individual curve
        # result.append('Summary of individual curves')
        # for k in range(len(self.curves)):
        #     result.append(str(k) + '/' + str(len(self.curves)))
        #     result.append(self.curves[k].summary())
        # result.append(separator)
        #
        # # summary of the cumulative statistics
        # result.append('Summary of the cummulative statistics')
        # result.append('pProt:\t\t' + str(self.pprot))
        # result.append('Contour length\tgamma\tks pValue')
        # for mu, gamma, pvalue in self.lengths:
        #     result.append(str(mu) + '\t\t' + str(gamma) + '\t' + str(pvalue))
        # result.append('Contour length histogram delimiting regions:')
        # result.append(str(self.hist_ranges))
        #
        # fname = self.info['name'] + '_results'
        # with open(fname, 'w') as myfile:
        #     myfile.write('\n'.join(result))
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
