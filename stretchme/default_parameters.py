import pandas as pd

default_parameters = {
    'p_prot': -1,
    'k_prot': -1,
    'l_prot': pd.DataFrame(),
    'p_dna': -1,
    'k_dna': -1,
    'l_dna': -1,
    'sheet_name': 0,
    'linker': None,
    'source': 'theory',
    'unit': 'nm',
    'speed': 500,    #nm/s
    'residues_distance': 0.365,
    'plot_columns': 2,
    'separator': ',',
    'states': None,
    'low_force_cutoff': 0.1,
    'significance': 0.0005,
    'max_distance': 0.3,
    'intervals': 1001,
    'spring_constant': 0.3,
    'bandwidth': 0.5,
    'max_force': 1600,
    'init_means': None,
    'method': 'marko-siggia',
    'initial_guess': {None: {'p_prot': 0.06, 'p_dna': 0.2, 'k_prot': 0.011, 'k_dna': 0.009, 'l_dna': 345, 'p_tot': 0.3,
                             'k_tot': 0.009, 'l_tot': 450},
                      'theory': {'p_prot': 0.3, 'p_dna': 0.13, 'k_prot': 0.008, 'k_dna': 0.009, 'l_dna': 345,
                                 'p_tot': 0.3, 'k_tot': 0.009, 'l_tot': 450},
                      'experiment': {'p_prot': 5.88, 'p_dna': 0.13, 'k_prot': 0, 'k_dna': 0.005, 'l_dna': 335,
                                     'p_tot': 0.3, 'k_tot': 0.009, 'l_tot': 450}}
 }
