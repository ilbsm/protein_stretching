from stretchme import Structure
import glob

directory = './'

experiment = Structure(debug=True, name='trmd_ca_4')
for file in glob.glob(directory + 'ca*.fd'):
    experiment.add_trace(file, columns=['D(1,N)', 'FORCE'], separator=' ', source='theory', linker='none', unit='A')
experiment.set_states(4)
experiment.analyze()
experiment.save_data()
