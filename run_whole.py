from stretchme import Structure
from simulation_data import data_partial
from os import remove, replace
import sys
import numpy as np

directory = '/Users/pawel/PycharmProjects/Rozciaganie/data/'


def run_whole(name, type, residues, missing, init_means=None, method='stretch-adjusted', state='unknotted'):
    name_type = name + '_' + type
    experiment = Structure(debug=True, name=name_type, state=state)
    file_data = data_partial[name_type]
    states_list = []
    run_parameters = {}
    for t in file_data:
        fname = directory + t['path']
        columns = t.get('columns', [])
        p_tot = t.get('p_prot', 0)
        k_tot = t.get('k_prot', 0)
        # p_dna = t.get('p_dna', 0)
        l_tot = t.get('l_dna', 0)
        # k_dna = t.get('k_dna', 0)
        l_linker = t.get('l_linker, 0')
        sheet_name = t.get('sheet_name', 0)
        case = [t.get('case', 0)]
        states = t.get('states', 3)
        states_list.append(states)
        if l_linker == 0:
            run_parameters = {
                'unit': 'A',
                'residues_distance': 0.38,
                'columns': columns,
                'separator': ' ',
                'p_prot': p_tot,
                'k_prot': k_tot,
                'method': method,
                'residues': residues,
                'cases': case,
                'missing': missing
            }
        else:
            run_parameters = {
                'residues_distance': 0.365,
                'sheet_name': sheet_name,
                'p_prot': 5.936,
                'k_prot': 0,
                # 'p_dna': p_dna,
                # 'k_dna': 0.009,
                # 'l_dna': 328.38,
                'method': method,
                'residues': residues,
                'missing': missing,
                'cases': case,
                'states': states,
                'linker': 'dna',
                'unit': 'nm',
                'state': state
            }
        experiment.add_trace(fname, **run_parameters)
        print(len(experiment.traces))
        # if len(experiment.traces) == 14:
        #      break
    experiment.set(overide=False, residues_distance=run_parameters['residues_distance'],
                   method=run_parameters['method'], residues=residues, num_states=max(states_list), init_means=init_means,
                   state=state)
    experiment.analyze()
    # experiment.save_data()
    return

def run_whole_experiment(name, method='stretch-adjusted'):
    name_type = name + '_exp'
    experiment = Structure(debug=True, name=name_type)
    file_data = data_partial[name_type]
    states_list = []
    run_parameters = {}
    for t in file_data:
        fname = directory + t['path']
        states = t.get('states', 3)
        states_list.append(states)
        run_parameters = {
            'residues_distance': 0.365,
            'sheet_name': t.get('sheet_name', 0),
            'p_tot': t.get('p_tot', 0),
            'k_tot': t.get('k_tot', 0),
            'l_tot': t.get('l_tot', 0),
            'l_linker': t.get('l_linker', 0),
            'first_range': t.get('first_range', (0,np.inf)),
            'last_range': t.get('last_range', (0, np.inf)),
            'missing': 14,
            'method': method,
            'residues': residues,
            'cases': [t.get('case', 0)],
            'states': states,
            'linker': 'dna',
        }
        experiment.add_trace(fname, **run_parameters)
        print(len(experiment.traces))
        # if len(experiment.traces) == 10:
        #     break
    experiment.set(overide=False, residues_distance=run_parameters['residues_distance'],
                   method=run_parameters['method'], residues=residues, num_states=max(states_list))
    experiment.analyze()
    return


proteins = {'5wyr': 248, 'trmd-no-knot': 240, 'trmd': 240, 'tm1570': 193, 'fuzja': 432}
missing = {'5wyr': 14, 'trmd-no-knot': 14, 'trmd': 14, 'tm1570': 14, 'fuzja': 28}
means = {('trmd', 'exp'): [45, 52, 87],
         ('trmd', 'kexp'): [40, 48, 82],
         ('tm1570', 'exp'): [10, 20, 25, 30, 35, 40, 45, 52, 57, 62, 69],
         ('tm1570', 'kexp'): [10, 20, 25, 30, 35, 40, 45, 52, 57, 62, 69],
         ('fuzja', 'exp'): None,
         ('fuzja', 'kexp'): None}

statess = {'exp': 'unknotted', 'kexp': 'knotted', 'pexp': 'partially'}
# protein, model = sys.argv[1], sys.argv[2]
to_c = [('exp', 'trmd'), ('kexp', 'tm1570'), ('exp', 'tm1570'), ('kexp', 'fuzja'), ('exp', 'fuzja'), ('pexp', 'fuzja'),
        ('cc', 'trmd'), ('cd', 'trmd'), ('ce', 'trmd'), ('pa', 'trmd'), ('cb', 'trmd'), ('pc', 'trmd'), ('pd', 'trmd'),
        ('sa', 'trmd'), ('sb', 'trmd'), ('sc', 'trmd'), ('sd', 'trmd'), ('ca', 'trmd-no-knot'), ('cb', 'trmd-no-knot'),
        ('cc', 'trmd-no-knot'), ('pa', 'trmd-no-knot'), ('pb', 'trmd-no-knot'), ('pc', 'trmd-no-knot'),
        ('sa', 'trmd-no-knot'), ('sb', 'trmd-no-knot'), ('sc', 'trmd-no-knot'), ('ca', '5wyr'), ('cb', '5wyr'),
        ('cc', '5wyr'), ('pa', '5wyr'), ('pb', '5wyr'), ('pc', '5wyr'), ('sa', '5wyr'), ('sb', '5wyr'), ('sc', '5wyr')]
# for model, protein in [('exp', 'trmd'), ('kexp', 'tm1570'), ('exp', 'tm1570'), ('kexp', 'fuzja'), ('exp', 'fuzja'), ('pexp', 'fuzja')]:
protein = 'trmd'
for model in ['exp']:
    print(model, protein)
    residues = proteins[protein]
    name_type = protein + '_' + model
    data_list = data_partial[name_type]
    run_whole_experiment(protein)


# for model in ['aa']:
#     print(model)
#     name_type = 'trmd' + '_' + model
#     data_list = data_partial[name_type]
#     run_whole('trmd', model)
#    replace('latex.info', model + '_whole_latex.info')
#    replace('data.info', model + '_whole_data.info')
