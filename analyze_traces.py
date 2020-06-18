from structure import Structure


''' 
The code analyzes the results from protein stretching. Each experiment, containing many stretching curves is expected to
be found in one file. The file list is given as the first parameter. The consecutive measurements are expected to be
found as ';' separated file, with consecutive columns d1;F1;d2;F2;d3;F3 etc...
The script calculates:
1. The rupture points;
2. Contour lengths;
3. The contour lengths histograms;
4. The total stretching;
5. Translates the [F,x] -> [F,L]
6. Histograms of contour lengths with numbers and positions of intermediate states;
7. Translates the [F,x] -> [x,L] with analysis of intermediate states;
8. Energy of transitions as the mean integral;
9. The histogram of rupture forces;
10. Lifetime of states from Dudko equation;
11. Parameters from Dudko equation;
12. Plots overlayed traces.

p.dabrowski@cent.uw.edu.pl
23.05.2020
'''

''' Parameters '''
data_files = [
    {'name': 'TrmD', 'residues': 240, 'distance': 4.4743, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm', 'speed': 1},
    # {'name': 'Tm1570', 'residues': 193, 'distance': 0.9901, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm', 'speed': 1},
    # {'name': 'fuzja', 'residues': 432, 'distance': 6.0558, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm', 'speed': 1},
    # {'name': 'trmd_CieplakT04', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_CieplakT05', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_CieplakT06', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_CieplakT05_spring', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_smogT04', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_smogT05', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_smogT06', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_smogT05_spring', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # {'name': 'trmd_CaT04', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_CaT05', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'trmd_CaT06', 'residues': 240, 'distance': 44.743, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CieplakT04', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CieplakT05', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CieplakT06', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CieplakT05_spring', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_smogT04', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_smogT05', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_smogT06', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_smogT05_spring', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CaT04', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CaT05', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': '5wyr_CaT06', 'residues': 248, 'distance': 51.034, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001},
    # {'name': 'AA_trmd', 'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory', 'speed': 0.001}
]
parameters = {
'data_path': '/Users/pawel/Documents/Projekty/2020 - Rozciaganie/data/',   # directory with the data
'data_file_prefix': 'raw_data_',                                      # the prefix for the core of the file name
'data_file_suffix': '.csv',                                           # the suffix for the core of the file name
'residues_distance': {'experiment': 3.65, 'theory': 3.86},            # distance between residues in stretched chain in A
'minimal_stretch_distance': 10,                                       # minimal distance between jumps
'high_force_cutoff': {'theory': 1, 'experiment': 5},                  # the cutoff delimiting the high force regime
'low_force_cutoff': {'theory': 0.1, 'experiment': 5},
'max_force_rupture': 42,                                              # the maximal force at which rupture may happen
'cluster_max_gap': 15,
'columns': 4,                                                         # number of columns of plots in mutliplots
'debug': True,
'initial_guess': {'pprot': 5.88, 'pdna': 0.16, 'K': 310, 'ldna': 345}
}
''' End of parameters '''


''' Main part '''
if __name__ == "__main__":
    for info in data_files:
        experiment = Structure(info, parameters, cases=[0])
        experiment.analyze()
