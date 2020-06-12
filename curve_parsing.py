""" The module containing tools for the analysis of a separate curve """

from tools import *
from datetime import datetime
from matplotlib import pyplot as plt


class Curve:
    def __init__(self, dist, forces, info, parameters, debug=False):
        self.dist = dist
        self.forces = forces
        self.info = info
        self.parameters = parameters
        self.debug = debug

        # smooth traces
        self.dist_smooth = []
        self.forces_smooth = []

        # derivatives of the traces
        self.dist_diff = []
        self.forces_diff = []

        # the ranges
        self.ranges = []

        # fitting
        self.coefficients = {}

        # rupture_forces
        self.rupture_forces = {}
        return

    def analyze(self):
        self._generate_smooths()
        self._generate_derivative()
        self._generate_ranges()
        self.fit()
        self.find_rupture_forces()
        self.find_energies()
        return

    # generating data
    def _generate_smooths(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Generating smooth traces.")
        errors = []
        f_smooth = []
        d_smooth = list(range(int(min(self.dist)), int(max(self.dist)) + 1))

        # averaging over the unit inteval
        for d in d_smooth:
            forces_interval = [self.forces[_] for _ in range(len(self.dist)) if d <= self.dist[_] < d + 1]
            try:
                f_smooth.append(sum(forces_interval) / len(forces_interval))
            except ZeroDivisionError:
                errors.append(d)

        # for the ones with nothing to average, corresponding distance is removed
        for d in list(reversed(errors)):
            d_smooth.pop(d_smooth.index(d))
        d_smooth, f_smooth = running_average(d_smooth, f_smooth)
        self.dist_smooth = d_smooth
        self.forces_smooth = f_smooth
        return

    def _generate_derivative(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Generating derivative of the traces.")
        d_diff, f_diff = find_derivative(self.dist_smooth, self.forces_smooth)
        self.dist_diff = d_diff
        self.forces_diff = f_diff
        return

    def _generate_ranges(self):
        """ Find the ranges between the ruptures. """
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Finding ranges.")
        ranges = []
        current_range = [self.dist_diff[0]]
        neg_diff = -1
        force_cutoff = self.parameters['high_force_cutoff'][self.info['source']]
        for d in range(len(self.dist_diff)):
            if self.forces_diff[d] < 0:
                if neg_diff < 0:
                    neg_diff = d
                # if the window of negative derivative is at least equal 'minimal_stretch_distance' parameter,
                # we finish the current range
                if neg_diff > 0 and self.dist_diff[d] - self.dist_diff[neg_diff] >= \
                        self.parameters['minimal_stretch_distance'] and len(current_range) == 1 \
                        and self.forces_smooth[d] >= force_cutoff:
                    current_range.append(self.dist_diff[neg_diff])
                    ranges.append(current_range)
                    current_range = []
            # establishing new range beginning as the first point when derivative becomes again positive
            if self.forces_diff[d] > 0:
                if neg_diff > 0 and len(current_range) == 0:
                    current_range.append(self.dist_diff[d])
                neg_diff = -1
        if len(current_range) == 1:
            current_range.append(self.dist_diff[-1])
            ranges.append(current_range)
        self.ranges = ranges
        # if self.debug:
        #     print("\t\tResult: " + '; '.join([str(round(x, 3)) + '-' + str(round(y, 3)) for x, y in self.ranges]))
        return

    # fitting and finding
    def fit(self):
        """ Finding the coefficients of WLC fit """
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Fitting the WLC model.")

        dist_to_fit = [np.array([self.dist[_] for _ in range(len(self.dist))
                                 if current_range[0] <= self.dist[_] <= current_range[1] and
                                 self.forces[_] > self.parameters['low_force_cutoff']]) for current_range in
                       self.ranges]
        forces_to_fit = [np.array([self.forces[_] for _ in range(len(self.dist))
                                   if current_range[0] <= self.dist[_] <= current_range[1] and
                                   self.forces[_] > self.parameters['low_force_cutoff']]) for current_range in
                         self.ranges]
        empties = []
        for k in range(len(dist_to_fit)):
            if len(dist_to_fit[k]) == 0:
                empties.append(k)
        for k in list(reversed(empties)):
            dist_to_fit.pop(k)
            forces_to_fit.pop(k)

        coefficients = fit_part_wlc(dist_to_fit[-1], forces_to_fit[-1])
        self._general_fit(dist_to_fit, forces_to_fit, coefficients)
        return

    def _general_fit(self, dist_array, forces_array, coefficients):
        # the model for symfit
        ds = []
        fs = []
        lengths = []
        pprot = Parameter('pprot', value=coefficients['pprot'], min=0.001)  # , min=3, max=10)
        model_dict = {}

        for k in range(len(dist_array)):
            ds.append(Variable('d' + str(k)))
            fs.append(Variable('F' + str(k)))
            lengths.append(Parameter('Lp' + str(k), value=max(dist_array[k])+20, min=max(dist_array[k])+1,
                                     max=coefficients['length']))
            model_dict[fs[k]] = pprot * (0.25 / ((1 - ds[k] / lengths[k]) ** 2) - 0.25 + ds[k]/lengths[k])

        model = Model(model_dict)
        arguments = {}
        for k in range(len(dist_array)):
            arguments['d' + str(k)] = dist_array[k]
            arguments['F' + str(k)] = forces_array[k]

        fit = Fit(model, **arguments)
        fit_result = fit.execute()

        # extracting coefficients
        self.coefficients['pprot'] = fit_result.value(pprot)
        self.coefficients['lengths'] = [fit_result.value(L) for L in lengths]
        return

    def find_rupture_forces(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Finding rupture forces.")
        # for r in range(len(self.ranges)):
        #     current_range = self.ranges[r]
        #     length = self.coefficients['L'][r]
        #     try:
        #         force = max([self.forces[_] for _ in range(len(self.dist))
        #                  if current_range[0] <= self.dist[_] <= current_range[1]])
        #     except ValueError:
        #         force = -1
        #     self.rupture_forces[length] = force
        return

    def find_energies(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Finding energies.")
        return

    # plotting
    def plot(self, position=None):
        position.set_xlim(min(self.dist), max(self.dist))
        position.set_ylim(0, max(self.forces))
        position.set_xlabel('Extension [A]')
        position.set_ylabel('Force [pN]')
        position.plot(np.array(self.dist), np.array(self.forces), label='Original data')
        position.plot(np.array(self.dist_smooth), np.array(self.forces_smooth), label='Smoothed data')

        pprot = self.coefficients['pprot']

        for r in range(len(self.ranges)):
            x0, x1 = self.ranges[r]
            position.hlines(max(self.forces) / 2, x0 + 2, x1 - 2, color=colors[r % len(colors)], alpha=0.5)
            position.axvline(x=x0, linestyle='--', color='#808080', linewidth=0.5)
            position.axvline(x=x1, linestyle='--', color='#808080', linewidth=0.5)

        for len_nr in range(len(self.coefficients['lengths'])):
            length = self.coefficients['lengths'][len_nr]
            d_plot = np.linspace(min(self.dist), length-1)
            f_plot = calc_forward_wlc(d_plot, length, pprot)
            label = 'Fit L=' + str(round(length, 3))   # + ' (' + str(residues) + ' AA)'

            position.plot(d_plot, f_plot, label=label, color=colors[len_nr % len(colors)])

        position.legend()
        return

    def summary(self):
        result = []
        result.append('Ranges:')
        result.append(str(self.ranges))
        result.append('pProt:\t\t' + str(self.coefficients['pprot']))
        result.append('Lengths:')
        result.append(str(self.coefficients['lengths']))
        return '\n'.join(result)
