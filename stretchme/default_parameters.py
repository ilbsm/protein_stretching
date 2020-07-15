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
    'speed': 1,
    'residues_distance': 0.365,
    'plot_columns': 4,
    'separator': ',',
    'states': None,
    'low_force_cutoff': 0.1,
    'significance': 0.001,
    'max_distance': 0.3,
    'intervals': 1001,
    'initial_guess': {None: {'p_prot': 0.7, 'p_dna': 0.1, 'k_prot': 0.005, 'k_dna': 0.005, 'l_dna': 345},
                      'theory': {'p_prot': 0.7, 'p_dna': 0.1, 'k_prot': 0.005, 'k_dna': 0.005, 'l_dna': 345},
                      'experiment': {'p_prot': 5.88, 'p_dna': 0.16, 'k_prot': 0, 'k_dna': 0.005, 'l_dna': 345}}
}