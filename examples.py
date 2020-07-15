from stretchme import Structure, simulate_experiment, analyze_experiment, analyze_trace
import glob

directory = ''


def analyze_file(filename, name, states=None):
    experiment = Structure(debug=True, name=name)
    for file in glob.glob(directory + filename):
        experiment.add_trace(file, columns=['D(1,N)', 'FORCE'], separator=' ', source='theory', linker=None, unit='A')
    experiment.set(states=states)
    experiment.analyze()
    experiment.save_data()


def simulate(traces=1, p_prot=0.7, k_prot=0.005, p_dna=0, k_dna=0, position_blur=0.1, force_blur=0.5, l_dna=350,
             l_prots=(25, 50, 100), rupture_forces=(5, 10), rupture_forces_blur=0.1,
             force_range=(0.1, 20), relaxation=0.1):
    result = simulate_experiment(traces=traces, p_prot=p_prot, k_prot=k_prot, p_dna=p_dna, k_dna=k_dna, l_prots=l_prots,
                                 position_blur=position_blur, force_blur=force_blur, force_range=force_range,
                                 relaxation=relaxation, rupture_forces=rupture_forces,
                                 rupture_forces_blur=rupture_forces_blur, l_dna=l_dna)
    result = result.sort_values(by=['d'])
    return result


def reverse_analyze(states=3):
    experiment = Structure(simulate(), source='theory', linker='none')
    experiment.set(states=states)
    experiment.analyze()
    return True


def trace_analysis(fname):
    analyze_trace(directory + fname, k_prot=0, p_prot=0.7, columns=['D(1,N)', 'FORCE'], separator=' ', unit='A',
                  debug=True, name='test' + fname)
    return True


trace_analysis('trmd_cg/ca06.fd')
