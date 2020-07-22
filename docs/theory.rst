.. _theory:

Theory
===========

.. _theory_fjc:

Freely-jointed chain
--------------------
One of the simplest models of the polymer is the freely-jointed chain.

.. _theory_wlc:

Worm-Like Chain model
---------------------
The work-like chain model may be viewed as the generalization of the freely-jointed chain model with the bond length :math:`b\rightarrow 0`.

.. _theory_fd:

The force-extension dependence in WLC model
-------------------------------------------
The stretching of the WLC polymer results in reduction of conformational space available for the polymer.

The exact form of the :math:`F(d)` dependency is not known. Most usually, the Marko-Siggia approximation is used (citation):

.. math::
    F(d;L,p) = p\left( \frac{1}{4\cdot(1-\frac{d}{L})^2} - \frac{1}{4} + \frac{d}{L} \right)
    :label: marko-siggia

The Marko-Siggia function commes from ...
Because of the construction, Marko-Siggia function agrees well with experiment in low- and high-force regime.
For the intermediate force up to 17% of discrepancy was reported (citation).

To decrease the error, some correction were proposed. In particular, the Petrosyan function:

.. math::
    F(d;L,p) = p\left( \frac{1}{4\cdot(1-\frac{d}{L})^2} - \frac{1}{4} + \frac{d}{L} - 0.8\cdot(\frac{d}{L})^{2.15} \right)
    :label: Petrosyan

was reported to decrease the error to 1% (citation).

.. _theory_ewlc:

Including the bond stretching
-----------------------------
The theory of the Worm-Like chain model assumes that the polymer is not extendible.
However, this is not the case in some real polymers. For example, upon stretching, the distances between the residues get stretched too.
The same is observed in case of *in silico* coarse-grained analysis, where the beads representing the polymer are confined with (extensible) harmonic spring.
To account for the stretching, the Elastic WLC model (eWLC) was introduced (citation).

In the eWLC, the expansion of the bond is included by introducing the elasticity :math:`k`.
Then, the relative extension :math:`x` is modified by the factor :math:`x \rightarrow x-\frac{F}{k}`:

.. math::
    F(x;p,k) = p\left( \frac{1}{4\cdot(1-x+\frac{F}{k})^2} - \frac{1}{4} + x - \frac{F}{k} \right)
    :label: reduced_ewlc

This simple modification is widely used e.g. in the analysis of stretching of DNA chains (citations).
However, such modification remains constant regardless of the chain contour length.
In fact, the chain with :math:`N_L` bonds, each of the same spring constant :math:`k`, under the force :math:`F` should exhibit additional expansion of :math:`N\cdot\frac{F}{k}`.
Then, the contour length (the maximal extension of the chain) increases as:

.. math::
    L \rightarrow L + N_L\cdot\frac{F}{k}
    :label: chain_length_increase

This should result in changing the WLC response as:

.. math::
    F(d;L,p,k) = p\left( \frac{1}{4\cdot(1-\frac{d}{L + N_L\cdot\frac{F}{k}})^2} - \frac{1}{4} + \frac{d}{L + N_L\cdot\frac{F}{k}} \right)
    :label: stretch_adjust1

This equation is rather unpleasant in manipulation, therefore, the extension can be expanded into series:

.. math::
    \frac{d}{L + N_L\cdot\frac{F}{k}} = \frac{d}{L} \frac{1}{1 + \frac{N_LF}{Lk}} \simeq x\cdot(1-\frac{N_L}{L}\frac{F}{k})
    :label: stretch_adjust_series

as from definition :math:`L = N_l \cdot b` where :math:`b` is the bond length, the equation :eq:`stretch_adjust1` can be presented as:

.. math::
    F(x;p,k) = p\left( \frac{1}{4\cdot(1-x\cdot\left(1-b\frac{F}{k})\right)^2} - \frac{1}{4} + x\cdot\left(1-b\frac{F}{k}\right) \right)
    :label: stretch_adjust2

This is the form of the **stretch-adjusted** wlc model introduced in the StretchMe package.

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
It is often needed to calculate the extension of the chain, knowing the force acting on it.
This problem also does not have a analytic solution. Usually, one neglects the behavior of the chain for small forces, and inverts the strong-force regime part, where the extension scales as square root of the force:

.. math::
    d(F;L,p_l) = L * \left(1-\frac{1}{2}\sqrt{\frac{Fp_L}{k_BT}}\right)
    :label: simple_inverting_marko_siggia

This approach is invalid outside the strong-force regime, e.g. when one protein state is ruptured in relatively low forces.
However, the *Marko-Siggia* function, after reduction to a polynomial is cubic in :math:`d` and therefore may be computationally easy inverted.
Namely, the value of :math:`d(F;L,p,k)` is the real root of the appropriate polynomial, lying in the interval :math:`(0,L)`.
The polynomials used are written explicitly in :ref:`theory_stretchme_functions`.

This approach is valid, as one can show, that the derivative of the *Marko-Siggia* function is strictly positive in the interval :math:`(0,L)`, mapping the interval into positive real numbers (for positive :math:`p` and :math:`k`).
Therefore, it is a bijection of these intervals, and therefore for each force there is only one :math:`d` within this interval.
In case of eWLC model, with non-zero elasticity :math:`k` one needs to use the implicit function theorem :ref:`theory_appendices_ift`.

.. _theory_contour_length:

Finding contour length
----------------------

Usually, the contours are found by global fitting of the curves.

.. _theory_knotting:

The expected chain reduction upon knotting
------------------------------------------

Upon knotting, the maximum length of the chain (contour length) gets decreased compared to the unknotted, stretched chain, as a portion of the chain is used to create a knot.

.. _theory_stretchme_functions:

The exact form of functions used in StretchMe
---------------------------------------------

In the computations, the form of the :math:`F(d)` dependencies are not suitable for calculations.
In particular, usually the basic quantity is the relative extension :math:`x=\frac{d}{L}`, but not the distance :math:`d`, nor the contour length :math:`L`.
Similarly, both in experiment and in the simulations, the temperature :math:`T` is constant, as well as the persistence length :math:`p_L`.
As a result, their quotient is constant and can be written as a single parameter :math:`p=\frac{k_BT}{p_L}`.
Finally, as the inextensible chain can be viewed as a chain can be viewed as an extensible chain with its parameter :math:`k\rightarrow\infty`.
This is, however, computationally uncomfortable and it is much better to operate with inverse of elasticity.

Therefore, the exact form of the methods used in the StretchMe package are as follows:

*Marko-Siggia*

.. math::
    F(x;p,k) = p\left( \frac{1}{4\cdot\left(1 - x + kF)\right)^2} - \frac{1}{4} + x - kF \right)
    :label: real_marko-siggia

This implicit function is calculated by forming a polynomial, cubic both in :math:`F` and in :math:`x`:

.. math::
    -\frac{F^3}{k^3}-\frac{F^3}{k^2 p}+\frac{3 F^2 x}{k^2}-\frac{2.25 F^2}{k^2}+\frac{2 F^2 x}{k p}-\frac{2 F^2}{k p}-\frac{3 F x^2}{k}+\frac{4.5 F x}{k}-\frac{1.5 F}{k}+\\
    \frac{F x^2}{p}+\frac{2 F x}{p}-\frac{F}{p}+x^3-2.25 x^2+1.5 x = 0
    :label: marko-siggia_polynomial


this allows for calculating both :math:`F(x;p,k)` as the sole real, positive root of this polynomial, as well as :math:`x(F;p,k)` as the sole real root lying within the interval :math:`(0,1)`.

*Stretch-adjusted*

.. math::
    F(x;p,k) = p\left( \frac{1}{4\cdot\left(1 - x\cdot(1-kF)\right)^2} - \frac{1}{4} + x\cdot(1-kF) \right)
    :label: real_stretch-adjusted

with the exact values of :math:`F(x;p,k)` and :math:`x(F;p,k)` calculated as the roots of the polynomial:

.. math::
    \frac{4 F^3 k^3 x^3}{b^3}+\frac{4 F^3 k^2 x^2}{b^2 p}-\frac{12 F^2 k^2 x^3}{b^2}+\frac{9 F^2 k^2 x^2}{b^2}-\frac{8 F^2 k x^2}{b p}+\frac{8 F^2 k x}{b p}+\frac{12 F k x^3}{b} + \\
    -\frac{18 F k x^2}{b}+\frac{6 F k x}{b}+\frac{4 F x^2}{p}-\frac{8 F x}{p}+\frac{4 F}{p}-4 x^3+9 x^2-6 x = 0
    :label: stretch-adjusted_polynomial

.. _theory_appendices:

Appendices
----------

.. _theory_appendices_ift:

Implicit Function Theorem
+++++++++++++++++++++++++

