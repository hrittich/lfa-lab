########
Glossary
########

.. vim: set spell spelllang=en_us:

Equations
=========

.. _poisson_equation:

Poisson Equation
----------------

Let the operator :math:`\Delta` be given by

.. math::

  \Delta u =
    \frac{\partial^2 u}{\partial x_1^2} + \cdots
    + \frac{\partial^2 u}{\partial x_d^2}
    \,,

where :math:`u` is a sufficiently smooth function.
The *Poisson equation* is

.. math:: -\Delta u = f \,,

where :math:`f` is a function.

Injection Operators
===================

.. _injection_restriction:

Injection Restriction
---------------------

Assume we have a fine grid :math:`\mathcal{G}_\mathbf{h}` and a coarse grid
:math:`\mathcal{G}_{\mathbf{h}'}` which is a subset of the fine grid.
The *injection restriction operator*
:math:`R_\mathrm{inj}: \ell_2(\mathcal{G}_\mathbf{h}) \to \ell_2(\mathcal{G}_{\mathbf{h}'})`
is defined by

.. math::

  [R_\mathrm{inj} u](x) = u(x)
  \quad \text{for all} \quad
  x \in \mathcal{G}_{\mathbf{h}'}
  \,.

.. _injection_interpolation:

Injection Interpolation
-----------------------

Assume we have a fine grid :math:`\mathcal{G}_\mathbf{h}` and a coarse grid
:math:`\mathcal{G}_{\mathbf{h}'}` which is a subset of the fine grid.
The *injection interpolation operator*
:math:`P_\mathrm{inj}: \ell_2(\mathcal{G}_{\mathbf{h}'}) \to \ell_2(\mathcal{G}_\mathbf{h})`
is defined by

.. math::

  [P_\mathrm{inj} u](x) =
  \begin{cases}
    u(x) & \text{for } x \in \mathcal{G}_{\mathbf{h}'} \\
    0    & \text{otherwise}
  \end{cases}
  \quad \text{for all} \quad
  x \in \mathcal{G}_{\mathbf{h}}
  \,.

Error Propagation Operators
===========================

.. _error_coarse_grid_correction:

Coarse Grid Correction
----------------------

The error propagator of the *coarse grid correction* is defined by

.. math::

   E = I - P (I - E_c) L_c^{-1} R L
   \,,

where

- :math:`L` is the linear system operator on the fine grid,
- :math:`L_c` is the linear system operator on the coarse grid,
- :math:`P` is the interpolation operator,
- :math:`R` is the restriction operator
- :math:`E_c` is the error propagator of the method that solves the coarse
  grid equation. In case of a two-grid method, :math:`E_c = 0`.

Software
========

.. _matplotlib:

Matplotlib
----------

Matplotlib is a Python library for visualizing mathematical functions and
data.

- Homepage: https://matplotlib.org/

