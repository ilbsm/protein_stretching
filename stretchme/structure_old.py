    def save_data(self):
        result = []
        separator = '################\n'

        # general info
        result.append('General info:')
        result.append('Name:\t\t\t' + self.info['name'])
        result.append('Residues:\t\t\t' + str(self.info['residues']))
        result.append('End-to-end distance:\t' + str(self.info['distance']))
        result.append('Linker:\t\t\t' + str(self.info['linker']))
        result.append('Unit:\t\t\t' + str(self.info['unit']))
        result.append('Data source:\t\t' + str(self.info['source']))
        result.append('Pulling speed:\t\t\t' + str(self.info['speed']))
        result.append(separator)

        # parameters
        result.append('Calculation parameters:')
        result.append('Data path:\t\t' + str(self.parameters['data_path']))
        result.append('Data file prefix:\t\t' + str(self.parameters['data_file_prefix']))
        result.append('Data file suffix:\t\t' + str(self.parameters['data_file_suffix']))
        result.append('Residue-residue distance:\t' + str(self.parameters['residues_distance']))
        result.append('Minimal distance between jumps:\t\t' + str(self.parameters['minimal_stretch_distance']))
        result.append('Low force cutoff:\t\t' + str(self.parameters['low_force_cutoff']))
        result.append('High force cutoff:\t\t' + str(self.parameters['high_force_cutoff']))
        result.append('Minimal gap between peaks in cluster:\t\t' + str(self.parameters['cluster_max_gap']))
        result.append('Number of columns in individual plots:\t\t' + str(self.parameters['columns']))
        result.append(separator)

        # summary of individual curve
        result.append('Summary of individual curves')
        for k in range(len(self.curves)):
            result.append(str(k) + '/' + str(len(self.curves)))
            result.append(self.curves[k].summary())
        result.append(separator)

        # summary of the cumulative statistics
        result.append('Summary of the cummulative statistics')
        result.append('pProt:\t\t' + str(self.pprot))
        result.append('Contour length\tgamma\tks pValue')
        for mu, gamma, pvalue in self.lengths:
            result.append(str(mu) + '\t\t' + str(gamma) + '\t' + str(pvalue))
        result.append('Contour length histogram delimiting regions:')
        result.append(str(self.hist_ranges))

        fname = self.info['name'] + '_results'
        with open(fname, 'w') as myfile:
            myfile.write('\n'.join(result))
        return

    def _calcualate_contour_lengths(self, pprot):
        self.contour_length = []
        for k in range(len(self.dists)):
            dist = np.array([self.dists[k][_] for _ in range(len(self.dists[k]))
                             if self.forces[k][_] > self.parameters['low_force_cutoff']])
            xs = np.array([invert_wlc_np(f, pprot) for f in self.forces[k] if
                           f > self.parameters['low_force_cutoff']])
            self.contour_length += list(dist / xs)
        self.contour_length = np.array(self.contour_length)
        return

    def _find_close_forces(self, r):
        forces = {}
        result = []
        current_range = self.hist_ranges[r]
        for c in self.curves:
            forces = merge_dicts(forces, c.rupture_forces)
        for force in forces.keys():
            if current_range[0] <= forces[force] <= current_range[1]:
                result.append(force)
        return result

    def _find_constants(self):
        pprots = [c.coefficients['pprot'] for c in self.curves]
        pprot = sum(pprots)/len(pprots)
        self._calcualate_contour_lengths(pprot)
        self.hist_ranges = find_hist_ranges(self.contour_length)
        for r in range(len(self.hist_ranges)):
            current_range = self.hist_ranges[r]
            values = [clength for clength in self.contour_length if current_range[0] <= clength <= current_range[1]]
            mu, gamma = cauchy.fit(values)
            test_statistic = cauchy.rvs(mu, gamma, len(values))
            self.lengths[mu] = {'gamma': gamma, 'pvalue':  ks_2samp(values, test_statistic).pvalue,
                                'forces': self._find_close_forces(r), 'ranges': current_range,
                                'color': colors[r % len(colors)]}
        self.minf = min([force for mu in self.lengths.keys() for force in self.lengths[mu]['forces']])
        self.maxf = max([force for mu in self.lengths.keys() for force in self.lengths[mu]['forces']])
        self.pprot = pprot
        return
