""" The stretchme package prepared for the analysis of the optical tweezers/AMF stretching of the biopolymers.
Creator: Pawel Dabrowski-Tumanski
e-mail: p.dabrowski-tumanski@uw.edu.pl
version 1.0"""

from .stretching_tools import *
from .structure import Structure
from .simulations import simulate_traces


def analyze_experiment(filename, experiment_name=None, sheet_name=0, separator=',', linker=None, unit='nm', speed=1.0,
                       source='theory', cases=None, read_columns=None, residues_distance=0.365, low_force_cutoff=0.1,
                       initial_guess=None, plot_columns=4, debug=False):
    """Main function for analysis all the data stored in one file. Creates the figures with each trace fitting, and the figure with analysis of all the data merged, including the histogram of rupture forces and the Dudko-Hummer-Szabo analysis of the states lifetimes.


        Args:
            filename (str/Pandas Dataframe): The path to the file with data, or the Pandas Dataframe with the data.
            experiment_name (str/None, optional): The name of the experiment, serving as the prefix for the files to be
                created. If None, the prefix will be created out of the input file, or it would be 'Experiment'.
                Default: None.
            sheet_name (str, optional): The name of the sheet in the .xls file to be analyzed. If 0, the first sheet
                will be analyzed. Default: 0.
            separator (str, optional): The separator delimiting the columns in the .csv file. Default: ','.
            linker (str/None, optional): The linker joining the protein to the beads which are pulled. The viable
                options are None or 'dna'. Default: None.
            unit (str, optional): The unit of the distance. Viable arguments are 'nm', or 'A'. Default: 'nm'
            speed (float, optional): The pulling speed, needed for Dudko-Hummer-Szabo analysis. Default: 1.0
            cases (list/None, optional): The list of cases to be analyzed from the data. If None, all the cases will be
                analyzed. Default: None.
            source (str, optional): The source of the data. Viable options are 'theory' or 'experiment'. The value used
                to establish the initial guess of the parameters.
            read_columns ([str, str]/None, optional): The headers of the columns in .xls/.csv files to be read. If None,
                every two columns will be treated as the distance and corresponding force. Default: None.
            residues_distance (float, optional): The expected distance in nm between Ca atoms of protein chain upon full
                stretching. Needed to calculate translate the contour lengths to the number of residues. Default: 0.365
            low_force_cutoff (float, optional): The lower cutoff of the force in pN to be analyzed. The lower the cutoff
                is set, the more outliers may appear. Default: 0.1.
            initial_guess (dict/None, optional): The dictionary with initial guesses for fitting the data. If None, the
                default initial_guess will be used. Default: None.
            plot_columns (int, optional): The number of traces in one row in the plots. Default: 4.
            debug (bool, optional): The debug mode. In the debug mode the log file is created. Default: False.

        Returns:
            bool: The return value. True for success, False otherwise.


        """

    parameters = pack_parameters(filename=filename, sheet_name=sheet_name, linker=linker, unit=unit, speed=speed,
                                 residues_distance=residues_distance, low_force_cutoff=low_force_cutoff, source=source,
                                 plot_columns=plot_columns, separator=separator, initial_guess=initial_guess)
    experiment = Structure(filename, cases=cases, columns=read_columns, parameters=parameters, name=experiment_name,
                           debug=debug)
    experiment.analyze()
    return True


def simulate_experiment(traces=1, p_prot=0.7, k_prot=0.005, l_prots=(25, 50, 100), p_dna=0, k_dna=0, l_dna=350,
                        rupture_forces=(10, 15), force_range=(0.1, 20), position_blur=0.1, force_blur=1,
                        rupture_forces_blur=0.1, relaxation=0.1):
    """Simulating the stretching data.


        Args:
            traces (int, optional): The number of traces to be simulated. Default: 1.
            p_prot (float, optional): The value kBT/pL for protein. Default: 0.7.
            k_prot (float, optional): The value of 1/K characterizing the protein. If 0, no effect of force on the
                distance between residues is seen. Default: 0.005.
            l_prots (tuple of floats, optional): The consecutive contour lengths of the protein during stretching.
                Default: (25, 50, 100).
            p_dna (float, optional): The value kBT/pL for DNA. If 0, no DNA is included. Default: 0.
            k_dna (float, optional): The value of 1/K characterizing the DNA. If 0, no effect of force on the
                distance between bases is seen. Default: 0.
            l_dna (float, optional): The contour length of the DNA. Default: 350.
            rupture_forces (tuple of floats, optional): The forces corresponding to the rupture of states determined by
                their contour length. The number of forces have to be number of states -1. Default: (10, 15).
            force_range ((float, float), optional): The range of forces to be measured. Default: (0.1, 20).
            position_blur (float, optional): The Gaussian width of the uncertainty of the position measurement.
                Default: 0.1.
            force_blur (float, optional): The Gaussian width of the uncertainty of the force measurement.
                Default: 1.
            rupture_forces_blur (float, optional): The Gaussian width of the uncertainty of the true rupture force.
                Default: 0.1.
            relaxation (float, optional): The parameter determining how long it takes the protein to relax after
                rupture, measured in the fraction of the contour length. Viable vales are between 0 and 1. Default: 0.1.

        Returns:
            Pandas Dataframe: Simulated distances and forces.


        """
    # TODO add a filter on the arguments

    return simulate_traces(traces=traces, p_prot=p_prot, k_prot=k_prot, l_prots=l_prots, p_dna=p_dna, k_dna=k_dna,
                           l_dna=l_dna, rupture_forces=rupture_forces, force_range=force_range,
                           position_blur=position_blur, force_blur=force_blur,  rupture_forces_blur=rupture_forces_blur,
                           relaxation=relaxation)


def analyze_trace(filename, case=0, columns=None, name=None, debug=False, **kwargs):
    """Analyzing single trace. Creates a figure with contour length histogram and fit.


        Args:
            filename (str/Pandas Dataframe): The path to the file with data, or the Pandas Dataframe with the data.
            case (int, optional): The case to be analyzed from the whole file. Default: 0
            columns: ([str, str]/None, optional): The headers of the distance and force columns in the file. If None,
                the first two columns will be used. Default: None.
            name: (str, optional): The name of the experiment, for the output data. If None, the name will be set to
                "Experiment". Default: None
            debug: (bool, optional): The debug mode. Default: False.
            **kwargs: other arguments passed to the reader or fitter. In particular, the fit parameters such as
                'p_prot', 'k_prot', 'p_dna', 'k_dna', 'l_dna' may be set.

        Returns:
            bool: True if successful, False otherwise.

        """
    experiment = Structure(filename, cases=[case], columns=columns, name=name, debug=debug, **kwargs)
    experiment.traces[0].analyze()
    experiment.traces[0].plot()
    return True
