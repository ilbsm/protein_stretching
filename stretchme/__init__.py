""" The stretchme package prepared for the analysis of the optical tweezers/AMF stretching of the biopolymers.
Creator: Pawel Dabrowski-Tumanski
e-mail: p.dabrowski-tumanski@uw.edu.pl
version 1.0"""

from .stretching_tools import *
from .structure import Structure
from .simulations import simulate_traces


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


def simulate_experiment(traces=1, p_prot=0.7, k_prot=200, p_dna=0, k_dna=None, position_blur=0.1, force_blur=1, l_dna=350,
                          l_prots=(25, 50, 100), rupture_forces=(10, 15), rupture_forces_blur=0.1,
                          force_range=(0.1, 20)):
    return simulate_traces(traces=traces, p_prot=p_prot, k_prot=k_prot, p_dna=p_dna, k_dna=k_dna, l_prots=l_prots,
                           position_blur=position_blur, force_blur=force_blur, force_range=force_range,
                           rupture_forces=rupture_forces, rupture_forces_blur=rupture_forces_blur, l_dna=l_dna)
