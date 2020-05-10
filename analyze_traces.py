from tools import extract_curve, find_ranges, fit_curve


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
10.05.2020
'''

''' Parameters '''
data_path = '/Users/pawel/Documents/Projekty/2020 - Rozciaganie/'   # directory with the data
data_files = {'TrmD': (240, 4.473, '3_1'),                          # dictionary with cores of file names as keys and
#              'Tm1570': (0, 0 , '3_1'),                             # a tuple with (number of residues between the
#              'fuzja': (0, 0, '3_1#3_1')}                           # attachment of DNA, distance in nm, expected knot)
                }
data_file_prefix = 'raw_data_'                                      # the prefix for the core of the file name
data_file_suffix = '.csv'                                           # the suffix for the core of the file name
residues_distance = 0.365                                           # distance between residues in stretched chain
gap_size = 0.2                                                      # the minimal gap between distances during jump
minimal_stretch_distance = 15                                       # minimal distance between jumps
break_size = 3                                                      # the minimal distance between the fitted parts
high_force_cutoff = 5                                               # the cutoff delimiting the high force regime

epsilon = 0.001
K_ranges = (250, 450)
''' End of parameters '''


def analyze_case(structure_name):
    coefficients = []
    for dist, forces in extract_curve(structure_name, data_path, data_file_prefix, data_file_suffix):
        ranges = find_ranges(dist, gap_size, minimal_stretch_distance, break_size)
        coefficients.append(fit_curve(dist, forces, ranges, high_force_cutoff))
        print(coefficients)
        expected = data_files[structure_name][0]*residues_distance
        missing = expected - coefficients['characteristic lengths (L)'][-1]
        missing_residues = missing/residues_distance
        print(expected, missing, missing_residues)
        break

    # for d, F in extract_curve(name):
    #     ranges = find_ranges(d)
    #     coefficients.append(fit_curve(d, F, ranges))
    #     Lx.append(translate_F(d, F, coefficients[-1]))
    #     find_energies(d, F, coefficients[-1])
    #     Fs += find_rupture_forces(d, F)
    # make_contour_length_histo(coefficients)
    # total_stretching = find_total_stretching(coefficients)
    # make_histograms(Lx)
    # make_force_histogram(Fs)
    # xd, G, v, tau = extract_life_times(Fs)
    return


''' Main part '''
for name in data_files.keys():
    analyze_case(name)
