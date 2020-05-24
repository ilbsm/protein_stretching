from tools import *


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
data_path = '/Users/pawel/Documents/Projekty/2020 - Rozciaganie/data/'   # directory with the data
# dictionary with cores of file names as keys and a tuple with (number of residues between the attachment of DNA, distance in nm, number of residues in linker between domains)
data_files = {
    # 'TrmD': {'residues': 240, 'distance': 4.4743, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm'},
    # 'Tm1570': {'residues': 193, 'distance': 0.9901, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm'},
    # 'fuzja': {'residues': 432, 'distance': 6.0558, 'linker': 'dna', 'source': 'experiment', 'unit': 'nm'},
    # 'CieplakT04': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CieplakT05': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CieplakT06': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CieplakT05_spring': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'smogT04': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'smogT05': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'smogT06': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    'smogT05_spring': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CaT04': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CaT05': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'},
    # 'CaT06': {'residues': 248, 'distance': 52.412, 'linker': 'none', 'unit': 'A', 'source': 'theory'}
}
data_file_prefix = 'raw_data_'                                      # the prefix for the core of the file name
data_file_suffix = '.csv'                                           # the suffix for the core of the file name
residues_distance = {'experiment': 3.65, 'theory': 3.86}            # distance between residues in stretched chain in A
gap_size = 0.2                                                      # the minimal gap between distances during jump
minimal_stretch_distance = 10                                       # minimal distance between jumps
break_size = 3                                                      # the minimal distance between the fitted parts
high_force_cutoff = 1                                               # the cutoff delimiting the high force regime
max_force_rupture = 42                                              # the maximal force at which rupture may happen
show_plots = True                                                   # boolean, if the plots are to be shown or saved
extension_speed = 200                                               # speed of extension in nm/s
force_decrease = 0.15                                               # decrease of smoothed force trace during rupture
high_force = 1                                                      # the force cut of the backgroud
columns = 4                                                         # number of columns of plots in mutliplots
''' End of parameters '''


def analyze_case(structure_name):
    coefficients = []
    rupture_forces = []
    k = 1
    contour_lengths = []
    ranges_total = []
    dist_total, forces_total = extract_curves(structure_name, data_path, data_file_prefix, data_file_suffix,
                                              unit=data_files[structure_name]['unit'])
    # dist_total_smooth = []
    # forces_total_smooth = []
    for dist, forces in zip(dist_total, forces_total):
        print(str(k) + '/' + str(len(dist_total)))
        ranges, dist_smooth, forces_smooth = find_ranges(dist, forces, gap=minimal_stretch_distance,
                                                         high_force_cutoff=high_force_cutoff)
        ranges_total.append(ranges)
        print(ranges)
        # ranges, dist_protein, forces_protein, coefficients_protein = separate_linker(ranges, dist_smooth,
        #                                                 forces_smooth, linker=data_files[structure_name]['linker'])
        if data_files[structure_name]['linker'] == 'dna':
            coefficients_trace = fit_curve_dna(ranges, dist_smooth, forces_smooth)
        else:
            coefficients_trace = fit_curve(ranges, dist_smooth, forces_smooth)
        coefficients.append(coefficients_trace)
        print(coefficients[-1])
        # rupture_forces += find_rupture_forces(dist, forces, ranges)
        # contour_lengths += transform_coordinates(dist, forces, coefficients[-1])
        # energies = find_energies(dist, forces, coefficients[-1])
        k += 1
    contour_length_gain = plot_coefficients(dist_total, forces_total, ranges_total, coefficients, name,
                           linker=data_files[structure_name]['linker'], show_plots=False, columns=columns,
                                            residues_distance=residues_distance[data_files[structure_name]['source']])

    # total_stretching = find_total_stretching(coefficients, show_plots)
    # force_counts, force_bins = make_force_histogram(rupture_forces, show_plots)
    # transformed_coordinates = plot_transformed_coordinates(contour_lengths, show_plots)
    # aver_dist, aver_forces = find_averages(dist_total, forces_total)
    # xd, g, v, tau = extract_life_times(force_counts, force_bins, aver_dist, aver_forces, extension_speed, show_plots)
    # save_data(contour_length_gain)
    return


''' Main part '''
for name in data_files.keys():
    analyze_case(name)
