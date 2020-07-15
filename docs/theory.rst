.. _theory:

Some theory
===========
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