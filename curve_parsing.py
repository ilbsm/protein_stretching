""" The module containing tools for the analysis of a separate curve """

from tools import *
from datetime import datetime


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
        self._generate_smooths()

        # derivatives of the traces
        self.dist_diff = []
        self.forces_diff = []
        self._generate_derivative()

        # the ranges
        self.ranges = []
        self._generate_ranges()

        # fitting
        self.coefficients = {}
        self.fit_errors = {}
        self.contourlengths = []

        # rupture_forces
        self.rupture_forces = {}
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
                # if the window of negative derivative is at least equal 'minimal_stretch_distance' parameter, we finish the current range
                if neg_diff > 0 and self.dist_diff[d] - self.dist_diff[neg_diff] >= self.parameters['minimal_stretch_distance'] \
                        and len(current_range) == 1 and self.forces_smooth[d] >= force_cutoff:
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
        if self.debug:
            print("\t\tResult: " + '; '.join([str(round(x, 3)) + '-' + str(round(y, 3)) for x, y in self.ranges]))
        return

    # fitting
    def fit(self):
        """ Finding the coefficients of WLC fit """
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Fitting the WLC model.")
        coefficients, errors = wlc_fit(self.ranges, self.dist_smooth, self.forces_smooth,
                                       self.info['linker'], self.parameters['low_force_cutoff'])
        self.coefficients = coefficients
        self.fit_errors = errors
        if self.debug:
            print("\t\tResult: " + str(self.coefficients))
            print("\t\tErrors: " + str(self.fit_errors))
        return

    def transform_coordinates(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Transforming coordinates to [F,L] space.")
        xs = np.array([invert_wlc_np(f, self.coefficients['pprot']) for f in self.forces if f > self.parameters['low_force_cutoff']])
        ds = np.array([self.dist[_] for _ in range(len(self.dist)) if self.forces[_] > self.parameters['low_force_cutoff']])
        if self.info['linker'] == 'dna':
            ddna = 2 * self.coefficients['ldna'] * np.array([invert_wlc_np(f, self.coefficients['pdna'])
                                                         for f in self.forces if f > 0])
            ds -= ddna
        self.contourlengths = list(ds/xs)
        return

    def find_rupture_forces(self):
        if self.debug:
            print("\t-> " + str(datetime.now().time()) + " Finding rupture forces.")
        for r in range(len(self.ranges)):
            current_range = self.ranges[r]
            length = self.coefficients['L'][r]
            try:
                force = max([self.forces[_] for _ in range(len(self.dist))
                         if current_range[0] <= self.dist[_] <= current_range[1]])
            except ValueError:
                force = -1
            self.rupture_forces[length] = force
        if self.debug:
            print("\t\tResult: " + str(self.rupture_forces))
        return

