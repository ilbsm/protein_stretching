default_parameters = {
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
    'initial_guess': {None: {'p_prot': 0.7, 'p_dna': 0.1, 'k_prot': 0.005, 'k_dna': 0.005, 'l_dna': 345},
                      'theory': {'p_prot': 0.7, 'p_dna': 0.1, 'k_prot': 0.005, 'k_dna': 0.005, 'l_dna': 345},
                      'experiment': {'p_prot': 5.88, 'p_dna': 0.16, 'k_prot': 0, 'k_dna': 0.005, 'l_dna': 345}}
}