########
Tutorial
########

For this introduction to LFA Lab, we assume that you have successfully
installed our software. (If this is not the case, see :ref:`installation`). We
start with a simple example.

.. _tut_jacobi_method:

The Jacobi Method
=================

For the beginning we start with the Jacobi method applied to the Poisson
equation (see :ref:`poisson_equation`). The following example uses LFA Lab to
compute spectral radius of the error propagator of the method.

.. literalinclude:: ../../demo/tutorial01-poisson.py

Let us walk through this example step by step.  To make any analysis, we first
have to define a grid, i.e., create an instance of
:py:class:`lfa_kernel.Grid`, because all operators that we can analyze by LFA
must be defined on a grid.

Then, we can create the operator corresponding to the Poisson equation. We
call the function :py:func:`lfa_lab.gallery.poisson_2d` for this purpose. The
module :py:mod:`lfa_lab.gallery` contains further predefined operators that
you can use.

In the next step, we create the error propagator of the Jacobi method for the
Poisson equation, using the :py:func:`lfa_lab.smoother.jacobi` function. This
function can only be used with operators that are defined by a stencil. How to
define such an operator will be discussed in the section
:ref:`defining_stencil_operators`.

The symbol of the error propagation operator of the Jacobi method is then
computed. Note, that this is the point where the actual computation happens,
all other operations before the computation of the symbol are merely
definitions. Since the computation of the symbol can be quite expensive, you
should store the symbol, when you want to investigate multiple of its
properties. (For the properties of a symbol, see
:py:class:`lfa_kernel.Symbol`.)

The last step is, to compute the spectral radius of the error propagation
operator. This command will output a value of 0.996 for the spectral radius.

Plotting Symbols
================

The *Fourier symbol* :math:`\hat{s}` of the operator :math:`S` is a
description of the operator :math:`S`.  The symbol assigns every value from
:math:`[0, \tfrac{2\pi}{h_1}) \times \cdots \times [0, \tfrac{2\pi}{h_d})` a
complex number, such that the following holds. If :math:`\hat{u}` and
:math:`\hat{f}` are the *discrete time Fourier transforms* of the functions
:math:`u` and :math:`f`, respectively, and :math:`f = S u` then

.. math:: \hat{f} = \hat{s} \cdot \hat{u}

(pointwise multiplication). Hence, by plotting the symbol of an operator we can
impove out understanding of the corresponding operator.

The following code plots the Fourer symbol of the error propagator of the
Jacobi method for the discrete Poisson equation.

.. literalinclude:: ../../demo/tutorial02-plot.py

This code is just a modification to the example from the previous section. To
plot the symbol, we import the :ref:`matplotlib` library. To plot the symbol
we use the function :py:func:`lfa_lab.plot.plot_2d`. (This function works only
for 2D problems. Problems in 1D should use the function
:py:func:`lfa_lab.plot.plot_2d`.) To actually display the plot, we need to use
the function `show()` of :ref:`matplotlib`. The result is shown below.

.. image:: tutorial02-plot.png

Note that LFA Lab uses the Frequency domain
:math:`[0, \tfrac{2\pi}{h_1}) \times \cdots \times [0, \tfrac{2\pi}{h_d})`.
This choice implies that the *low modes* are located near the corners of the
plot. (Recall that every mode that is node a low mode is called a *high
mode*).

By inspecting the plot and using the relation
:math:`\hat{f} = \hat{s} \cdot \hat{u}`, we see that applying the operator
:math:`S` returns a function that, in comparison to the input functions, has
smaller function values on the high modes. Hence, after some applications of
:math:`S` the low modes will be dominating.

It is known that a function whose discrete time Fourier transform is dominant
in the low modes is slowly varying. Thus, multiple applications of :math:`S`
result in a function which is slowly varying. Therefore, the Jacobi method is
a *smoother* for the discrete Poisson equation.

Smoothing Analysis
==================

The smoothing effect of an operator :math:`S` can be quantified by a so called
*smoothing analysis*. The idea of this analysis is to choose a filtering
operator :math:`Z` that maps the low modes to zero while leaving the high
modes unchanged. (An operator like that is also called an *idealized coarse
grid correction*.)
Then, :math:`ZS` characterizes how :math:`S` acts on the high
modes. Hence, the smoothing effect of the operator :math:`ZS` is quantified by
computing the spectral radius

.. math::

  r(ZA)
  \,.

The following code performs this computation.

.. literalinclude:: ../../demo/tutorial03-smoothing.py


.. _tutorial_two_grids:

Two Grids
=========

The two-grid method requires

- a smoothing operator :math:`S` and
- a coarse grid correction :math:`E_\mathrm{cgc}`.

The error propagator of the two-grid method is the given by

.. math:: E = S E_\mathrm{cgc} S \,.

To construct the coarse grid correction error propagator we need

- the linear system operator of the fine grid `L`,
- the linear system operator of the coarse grid `Lc`,
- the restriction operator `R`, and
- the interpolation operator `P`.

This error propagator can be computed using the
:py:func:`lfa_lab.two_grid.coarse_grid_correction` function. This function
does nothing more than to evaluate the formula (see
:ref:`error_coarse_grid_correction`) for the coarse grid correction. It is,
however, better to use the predefined function in LFA lab to avoid typos.

The following code computes the required operators, then the error propagator
of the coarse grid correction, and then the error propagator of the entire
two-grid method.

.. literalinclude:: ../../demo/tutorial04-twogrid.py


.. _defining_stencil_operators:

Defining Stencil Operators
==========================

Up to this point we used the functions from the
:py:mod:`lfa_lab.gallery` module, to create different operators. We shall now
discuss how to define these operators manually.
Let us start with the Laplace operator.

The stencil of the discrete Laplace operator is

.. math::

  \begin{bmatrix}
    & -\frac{1}{h_2^2} \\
    -\frac{1}{h_1^2} & \left( \frac{2}{h_1^2} + \frac{2}{h_2^2} \right) & -\frac{1}{h_1^2} \\
    & -\frac{1}{h_2^2}
  \end{bmatrix}_{\mathbf{h}}
  \,.

In LFA Lab stencils are described by a list containing stencil entries. Each
stencil entry is a tuple consisting of the offset of the entry and the
coefficient of the entry::

  [ ((o1_x, o1_y, ...), v1), ((o2_x, o2_y, ...), v2), ... ]

Using the :py:func:`lfa_lab.operator.from_stencil` function, we can turn this
description of a stencil into an operator that we can use in our analysis. The
following code constructs the stencil for the discrete Laplace operator::

  h1 = 1.0/32
  h2 = 1.0/32
  entries = [
      (( 0, -1), -1.0 / (h2*h2)),
      ((-1,  0), -1.0 / (h1*h1)),
      (( 0,  0),  2.0 / (h1*h1) + 2.0 / (h2*h2)),
      (( 1,  0), -1.0 / (h1*h1)),
      (( 0,  1), -1.0 / (h2*h2))
    ]
  L = operator.from_stencil(entries, grid)

Defining an interpolation and restriction operator is slightly more
complicated.

A stencil operator defined using :py:func:`lfa_lab.operator.from_stencil` maps
from a grid with a step-size :math:`\mathbf{h}` to a grid with the same
step-size :math:`\mathbf{h}`. Combining a stencil operator with a canonical
injection or canonical restriction (see :ref:`injection_restriction` and
:ref:`injection_interpolation`), however, is usually sufficient
to define the desired interpolation or restriction operator.

More precisely, to define an interpolation operator we usually combine a
stencil operator :math:`P_\mathrm{st}` and the injection interpolation
:math:`P_\mathrm{inj}` to define the desired interpolation by

.. math::

  P = P_\mathrm{st} P_\mathrm{inj}
  \,.

Furthermore, to define a restriction operator we usually combin a stencil
operator :math:`R_\mathrm{st}` and the injection restriction
:math:`R_\mathrm{inj}` to define the destired restriction by

.. math::

   R = R_\mathrm{inj} R_\mathrm{st}
   \,.

For example to construct the linear interpolation with standard coarsening in
2D we can use the following code::

  entries = [
      ((-1, -1), 1.0/4),
      (( 0, -1), 1.0/2),
      (( 1, -1), 1.0/4),
      ((-1,  0), 1.0/2),
      (( 0,  0), 1.0),
      (( 1,  0), 1.0/2),
      ((-1,  1), 1.0/4),
      (( 0,  1), 1.0/2),
      (( 1,  1), 1.0/4),
    ]
  P = operator.from_stencil(entries, fine) * \
          operator.injection_interpolation(fine, coarse)

We give below a complete example for user defined stencils.

.. literalinclude:: ../../demo/tutorial05-stencils.py


Multigrid Method
================

To define the error propagation operator of a multigrid method, we use the
`coarse_error` argument of the
:py:func:`lfa_lab.twogrid.coarse_grid_correction` function. This argument
gives the error propagator of the coarse grid solver. Hence, we can write a
recursive function that computes the error propagator of the multigrid method.

.. literalinclude:: ../../demo/tutorial06-multigrid.py


