.. _theory:

Theory
===========

.. _theory_fjc:

Freely-jointed chain
--------------------

.. _theory_wlc:

Worm-Like Chain model
---------------------

.. _theory_fd:

The force-extension dependence in WLC model
-------------------------------------------
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

This can be given in the reduced form, for a fixed :math:`L`, with :math:`x=d/L`:

.. math::
    F(x;p,k) = p\left( \frac{1}{4\cdot(1-x+kF)^2} - \frac{1}{4} + x - kF \right)
    :label: reduced_ewlc


When :math:`k=0` one recovers the original Marko-Siggia function.

.. _theory_ewlc:

Extensible WLC
--------------

.. _theory_linker:

Influence of the linker
-----------------------
It is common, that the polymer of interest (e.g. protein) is attached to the beads which are pulled with a linker.
The linker itself can also be a polymer (e.g. DNA), which could be modelled again by WLC model.
In such a case, the system may be compared to a series of springs.
As the force on each two ends of the spring is equal, the same force :math:`F` stretches each spring.
As a result, the total distension :math:`d` registered in the system of a protein connected by two DNA linkers is equal:

.. math::

    d(F) = d_{DNA_1}(F,L_{DNA_1},p_{DNA_1},{DNA_1}) + d_{PROT}(F,L_{PROT},p_{PROT},k_{PROT}) + d_{DNA_2}(F,l_{DNA_2},p_{DNA_2},k_{DNA_2})

where :math:`L` is the contour length, :math:`p` persistence length, and :math:`k` is elasticity.
Usually, the DNA linkers on each side are identical in structure, therefore, their distensions are equal and one can neglect the last term, introducing the 'effective' parameters for DNA:

.. math::

    d(F) = d_{DNA}(F,L_{DNA},p_{DNA},{DNA}) + d_{PROT}(F,L_{PROT;i},p_{PROT},k_{PROT}) + d_{DNA}(F,l_{DNA},p_{DNA},k_{DNA})


.. _theory_inverting:

Inverting the force-extension dependence
----------------------------------------
The need to include the effect of the linker requires inverting the force-extension dependence.
As in the equation :eq:`reduced_ewlc` the force is given in the implicit function, it cannot be simply inverted.
However, one can calculate its derivative from the :ref:`theory_appendices_ift`:

To see that it is positive in for :math:`0<x<1` and any :math:`p,k>0`.
Therefore, this function is reversible in this region.
By rearanging the terms, one can show, that the function :eq:`reduced_ewlc` is equivalent to the following:


.. _theory_contour_length:

Finding contour length
----------------------

.. _theory_knotting:

The expected chain reduction upon knotting
------------------------------------------

.. _theory_appendices:

Appendices
----------

.. _theory_appendices_ift:

Implicit Function Theorem
+++++++++++++++++++++++++

