from .stretching_tools import *
from .trace import Trace
from datetime import datetime


class Structure:
    def __init__(self, filename, cases, parameters, name=None, debug=False):
        if debug:
            print("Analyzing\t" + filename + "\nCases:\t" + str(cases) + "\nParameters:\t" + str(parameters))
        self.data = read_data(filename, cases, sheet_name=parameters['sheet_name'], separator=parameters['separator'])
        self.parameters = parameters
        self.debug = debug
        self.traces = []
        self.init_time = datetime.now().time()
        if name:
            self.name = name
        else:
            self.name = filename.split('/')[-1].split('.')[0]

        if self.debug:
            print(str(datetime.now().time()) + "\t Found " + str(int(len(list(self.data))/2)) + " traces to analyze.")
        for number in list(set([name.strip('d_F') for name in list(self.data)])):
            self.traces.append(Trace(number, self.data.loc[:, ['d_' + number, 'F_' + number]],
                                     parameters=self.parameters, debug=self.debug))
        return

    def _find_constants(self):
        return

    def _make_partial_plots(self):
        if self.debug:
            print(str(datetime.now().time()) + " Making plots of individual fits.")
        number = len(self.traces)       # number of traces to plot
        columns = self.parameters['columns']
        rows = max(int(np.ceil(float(2*number) / columns)), 2)
        fig, axes = plt.subplots(rows, columns, dpi=600, figsize=(5*int(columns), 5*int(rows)))
        max_contour_length = max([t.max_contour_length for t in self.traces])
        if max_contour_length <= 0:
            max_contour_length = None
        for k in range(0, number):
            self.traces[k].plot_histogram(position=axes[int(2*k/4), (2*k) % 4], max_contour_length=max_contour_length)
            self.traces[k].plot_fits(position=axes[int((2*k+1)/4), (2*k+1) % 4])
        fig.tight_layout()
        if self.debug:
            print("Saving contour lengths figure to " + self.name + '_contour_lengths.png')
        plt.savefig(self.name + '_contour_lengths.png')
        return

    def _plot_individual_contour_length_histo(self, position):
        position.set_xlabel('Fitted contour length')
        position.set_ylabel('Occurences')
        position.set_title('Histogram of fitted contour lengths')
        lengths = [parameter[0] for c in self.curves for parameter in c.coefficients['lengths']]
        position.hist(lengths, bins=50)
        return

    def _plot_total_contour_length_histo(self, position):
        # lengths_range = np.linspace(min(self.contour_length), max(self.contour_length))
        position.set_xlabel('Contour length')
        position.set_ylabel('Occurences')
        position.set_title('Contour length histogram')
        hist_values = []
        for t in self.traces:
            hist_values += list(t.data['hist_values'])
        position.hist(hist_values, bins=200, density=True)
        return

    def _plot_forces_histogram(self, position):
        position.set_xlabel('Force [pN]')
        position.set_ylabel('Occurences')
        position.set_title('Rupture force histogram')
        # position.set_xlim(self.minf, self.maxf)
        # fspace = np.linspace(self.minf, self.maxf, 100)
        # for mu in self.lengths.keys():
        #     forces = self.lengths[mu]['forces']
        #     n, bins, patches = position.hist(forces, bins=len(forces)/4, density=True,
        #                                      facecolor=self.lengths[mu]['color'], alpha=0.5)
        #     mu, sigma = norm.fit(forces)
        #     label = str(round(mu, 3))
        #     position.plot(fspace, norm.pdf(fspace, mu, sigma), self.lengths[mu]['color'], linestyle='--', label=label)
        #     self.lengths[mu]['hist_vals'] = n
        #     self.lengths[mu]['hist_width'] = bins[1] - bins[0]
        #     self.lengths[mu]['hist_beg'] = bins[0]
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
        if self.debug:
            print(str(datetime.now().time()) + " Making histograms.")
        fig, axes = plt.subplots(2, 2, dpi=600, figsize=(10, 10))

        # the histogram of contour plots from individual traces
        self._plot_individual_contour_length_histo(axes[0, 0])

        # the total contour length histogram
        self._plot_total_contour_length_histo(axes[0, 1])

        # the rupture forces histogram
        self._plot_forces_histogram(axes[1, 0])

        # Dudko analysis
        self._plot_dudko_analysis(axes[1, 1])

        fig.tight_layout()
        if self.debug:
            print("Saving histograms figure to " + self.name + '_histograms.png')
        plt.savefig(self.name + '_histograms.png')
        return

    def make_plots(self):
        self._make_partial_plots()
        self._make_histograms()
        return

    def save_data(self):
        return

    def add_trace(self, filename):
        return

    def analyze(self):
        for t in self.traces:
            t.analyze()
        self._find_constants()
        self.make_plots()
        self.save_data()
        return
