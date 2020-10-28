from stretchme import analyze_trace
from simulation_data import data_partial
from os import remove, replace
import sys

directory = '/Users/pawel/PycharmProjects/Rozciaganie/data/'


def run_case(name, type, case, method='stretch-adjusted', **kwargs):
    name_type = name + '_' + type
    file_data = data_partial[name_type][case]
    fname = directory + file_data['path']
    # states = file_data['states']
    columns = file_data.get('columns', [])
    sheet_name = file_data.get('sheet_name', 0)
    parameters = {**file_data, **kwargs}
    p_prot = parameters.get('p_prot', 0)
    k_prot = parameters.get('k_prot', 0)
    p_dna = parameters.get('p_dna', 0)
    k_dna = parameters.get('k_dna', 0)
    l_dna = parameters.get('l_dna', 0)
    first_range = parameters.get('first_range', None)
    last_range = parameters.get('last_range', None)
    missing = parameters.get('missing', 0)
    residues = parameters.get('residues', 240)
    state = parameters.get('state', 'unknotted')
    first_range = parameters.get('first_range', None)
    last_range = parameters.get('last_range', None)

    exp_name = 'results/'
    if l_dna > 0:
        exp_name += 'experiment/'
    exp_name += name + '/'
    if not state == 'knotted':
        exp_name += 'not_decomposed/'
    exp_name += '_'.join([name_type, str(case+1), str(missing)])


    if l_dna == 0:
        # exp_name = 'results/' + name + '/' + '_'.join([name_type, str(case), str(p_prot), str(k_prot)])
        run_parameters = {'residues_distance': 0.38,
                          'columns': columns,
                          # 'p_prot': p_prot,
                          # 'k_prot': k_prot,
                          'unit': 'A',
                          'separator': ' ',
        #                  'states': states,
                          'trace_name': '_'.join([name, type, str(case + 1)]),
                          'method': method,
                          'name': exp_name,
                          'residues': residues,
                          'missing': missing,
                          'state': state}
    else:
        # exp_name = 'results/' + name + '/' + '_'.join([name_type, str(case), str(p_prot), str(p_dna), str(k_dna),
        #                                                str(l_dna)])
        run_parameters = {'residues_distance': 0.365,
                          'sheet_name': sheet_name,
                          # 'p_prot': 5.936,
                          # 'k_prot': 0.009,
                          # 'p_dna': p_dna,
                          # 'k_dna': 0.009,
                          # 'l_dna': l_dna,
                          # 'states': states,
                          'method': method,
                          'trace_name': '_'.join([name, type, str(case + 1)]),
                          'case': [parameters.get('case', 0)],
                          'name': exp_name,
                          # 'linker': 'dna',
                          'missing': missing,
                          'residues': residues,
                          'state': state,
                          'first_range': first_range,
                          'last_range': last_range}

    analyze_trace(fname, debug=True, **run_parameters)
    return True

# for l_dna in [345]:
#     for p_prot in [6.0]:
#         for p_dna in [0.16]:
#             for k_dna in [0.0088]:
#                 print(p_prot, p_dna, k_dna)
#                 try:
#                     run_case('experiment', 'fuzja', 0, k_dna=k_dna, p_dna=p_dna, p_prot=p_prot, l_dna=l_dna)
#                 except ValueError:
#                     print("I've got problem")


proteins = {'5wyr': 248, 'trmd-no-knot': 240, 'trmd': 240, 'tm1570': 193, 'fuzja': 432}
# protein, model, case, missing = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
protein, model, case, missing = 'trmd', 'exp', 1, 14
if len(sys.argv) > 5 and int(sys.argv[5]) > 0:
    state = 'knotted'
else:
    state = 'unknotted'
# protein, model, case, missing = 'fuzja', 'exp', 17, 28
print(model, protein, case)
residues = proteins[protein]
name_type = protein + '_' + model
data_list = data_partial[name_type]
run_case(protein, model, case, missing=missing, residues=residues, state=state)

