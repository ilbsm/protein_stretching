""" The stretchme package prepared for the analysis of the optical tweezers/AMF stretching of the biopolymers.
Creator: Pawel Dabrowski-Tumanski
e-mail: p.dabrowski-tumanski@uw.edu.pl
version 1.0"""

from .stretching_tools import *
from .structure import Structure


def analyze_experiment(filename, sheet_name=0, residues=None, distance=None, linker=None, source=None, unit=None,
                       speed=None, residues_distance=0.365, minimal_stretch_distance=3, high_force_cutoff=None,
                       cases=None, low_force_cutoff=None, max_rupture_force=None, max_cluster_gap=15, plot_columns=4,
                       initial_guess=None, separator=',', output=None, read_columns=None, debug=True):
    parameters = pack_parameters(filename=filename, sheet_name=sheet_name, residues=residues, distance=distance,
                                 linker=linker, source=source, unit=unit, speed=speed,
                                 residues_distance=residues_distance, minimal_stretch_distance=minimal_stretch_distance,
                                 high_force_cutoff=high_force_cutoff, low_force_cutoff=low_force_cutoff,
                                 max_rupture_force=max_rupture_force, max_cluster_gap=max_cluster_gap,
                                 plot_columns=plot_columns, separator=separator, initial_guess=initial_guess)
    experiment = Structure(filename, cases=cases, columns=read_columns, parameters=parameters, name=output, debug=debug)
    experiment.analyze()
    return


def add_trace(structure, filename):
    return structure.add_trace(filename)
