.. _Overview:

Stretchme overview
==================
Stretchme is a Python3 package allowing for the automated analysis of the stretching of proteins with optical tweezers or AFM. The expected input from the user are the data *(d, F)* with *d* being the stretching distance and *F* the measured force.

Strechme analyses data both from simulations (theory) and experiment, with or without DNA linker joining the protein to the pulled bead. Strechme can also analyze the data, where the protein is being stretched simultaneously with the DNA.

Basically, as a result, the user should obtain two figures. One figure contains all the cumulative data for the experiment, including the histograms of obtained states fits to the overlayed traces, the histograms of forces, and the Dudko-Hummer-Szabo analysis.

.. figure::  images/example_whole.jpg
   :align:   center

   The exemplary output of the cumulative data.

The second figure contains the histogram of observed states, separated by their contour length, along with the fits for individual measurements (traces).

.. figure::  images/example_single.jpg
   :align:   center

   The exemplary output of the analysis of single trace.
