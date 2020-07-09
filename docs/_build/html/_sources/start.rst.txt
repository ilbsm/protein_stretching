.. _Start:

Quick start
===========

.. _start_installation:

Installation
------------
The installation of the package is as simple as writing::

    python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps strechme

It is distributed (yet) at the test.pypi.org::

    https://test.pypi.org/project/Stretchme

As this is still a test version, it is maintained at a test server. In nearest future it will be provided at the real pypi.
.. _start_theory:

Some theory
-----------
The force-extension curves for proteins and DNA are described using a version of elastic worm-like chain model. In the
eWLC model, the force *F* dependence on the extension *d* is given by implicit equation:

.. math::

    F(d;L,p_l,k_l,T) = \frac{k_BT}{p_l}\left( \frac{1}{4\cdot(1-\frac{d}{L}+\frac{F}{k_l})^2} - \frac{1}{4} + \frac{d}{L} - \frac{F}{k_l} \right)

which depends parametrically on the protein contour length (:math:`L`), persistence length (:math:`p_l`), stretching
constant of the bonds (:math:`k_l`) and the temperature (:math:`T`).
As usually the temperature and persistence length are constant during experiment, one can treat their quotient as
another variable :math:`p=\frac{k_BT}{p_l}`. Similarly, as the values of :math:`k_l` are usually high, in the everyday
application it is more convenient to work with its inverse :math:`k=\frac{1}{k_l}`. This yields the equation which is
fitted in the Strechme package:

.. math::

    F(d;L,p,k) = p\left( \frac{1}{4\cdot(1-\frac{d}{L}+kF)^2} - \frac{1}{4} + \frac{d}{L} - kF \right)

When :math:`k=0` one recovers the original Marko-Siggia function.

.. _start_requirements:

Requirements
------------
The package runs for Python 3. In requires the following packages, installed during the standard package installation:
* setuptools
* pandas
* matplotlib
* numpy
* scipy
* scikit-learn
* lmfit

.. _start_examples:

Example usage
-------------
Below the most common usages of the package are shown. For tweaking the analysis, see the documentation of the main classes *Structure* and *Trace*.

Example 1
+++++++++
Analyzing the 'data.xls' file with the measured distances and corresponding forces in consecutive columns is very easy::

    from stretchme import analyze_experiment
    analyze_experiment('data.xls')

results with the complete analysis of the data. One can also specify the sheet name::

    from stretchme import analyze_experiment
    analyze_experiment('data.xls', sheet_name='my_first_protein')

and tell Strechme, that the data come from *experiment* with *dna* linkers attaching the protein to the beads::

    from stretchme import analyze_experiment
    analyze_experiment('data.xls', sheet_name='my_first_protein', source='experiment', linker='dna')

Example 2
+++++++++
One can also analyze traces stored in specific files, by manipulating the parameters of the *Structure* class::

    from stretchme import Structure
    experiment = Structure()
    for file in glob.glob('my_proteins*.csv):
        experiment.add_trace(file)
    experiment.analyze()

The last command will start the whole analysis, resulting in the histograms and fitted traces.

In general, the data in the files are expected to contain distance and corresponding forces. If columns with specified headers (e.g. 'stretch_distance' and 'measured_force') are to be read, this can be specified within the *add_trace* method::

    from stretchme import analyze_experiment
    experiment = Structure()
    for file in glob.glob('my_proteins*.csv):
        experiment.add_trace(file, columns=['stretch_distance', 'measured_force'])
    experiment.analyze()

One can also specify the expected number of states (e.g. 4) in the whole experiment::

    from stretchme import analyze_experiment
    experiment = Structure()
    for file in glob.glob('my_proteins*.csv):
        experiment.add_trace(file, columns=['stretch_distance', 'measured_force'])
    experiment.set(states=4)
    experiment.analyze()

This will set the value of expected states in the total analysis for 4, but for each individual case, the number of states may still be different. If you want to set the individual trace, this can be specified by invoking the *.set()* method on the entries of the *Structure.traces* list::

    from stretchme import analyze_experiment
    experiment = Structure()
    for file in glob.glob('my_proteins*.csv):
        experiment.add_trace(file, columns=['stretch_distance', 'measured_force'])
        experiment.traces[-1].set(states=4)
    experiment.analyze()

Example 3
+++++++++
In general, fitting the curves takes most of the time. If you know the values of *p_prot*, *k_prot* etc, you can set
them before analysis, prohibiting therefore finding new parameters. You can set the parameters globally::

    experiment = Structure('data.xls')
    experiment.set(k_prot=0.005, p_prot=0.7)
    experiment.analyze()

 or for each trace individually::

    experiment = Structure()
    for file in glob.glob('my_proteins*.csv):
        experiment.add_trace(file, columns=['stretch_distance', 'measured_force'])
        experiment.traces[-1].set(p_prot=0.7, k_prot=0.005)
    experiment.analyze()

Example 4
+++++++++
One can also simulate the traces and plot them using standard MatPlotLib functions. The data obtained are as raw as in experiment, i.e. they are not sorted, therefore, in order to obtain a good plot, one needs to sort them first::

    from stretchme import simulate_experiment
    from matplotlib import pyplot as plt
    data = simulate_experiment(traces=1)
    data = data.sort_values(by=['d'])
    plt.plot(data['d'], data['F'])
    plt.show()

