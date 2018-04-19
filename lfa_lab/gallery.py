# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018  Hannah Rittich
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from __future__ import division

from .core import *
from .stencil import *
from .dag import *
from . import operator
from .util import *
import numpy as np


__all__ = [
    'fw_restriction',
    'ml_interpolation_stencil',
    'ml_interpolation',
    'fw_restriction_stencil',
    'poisson_2d',
]

def ml_interpolation_stencil(fine_grid, coarse_grid):
    """Return the stencil of a multilinear interpolation.

    :rtype: SparseStencil
    """
    d = fine_grid.dimension()
    ones = np.ones(d)
    zeros = np.zeros(d)

    coarsening = np.array(fine_grid.coarsening_factor(coarse_grid))

    positions = ShiftedGrid.from_corners(
            -coarsening+ones,
             coarsening-ones)

    # the 1d linear basis function
    def basis_1d(x, n):
        if x >= 0:
            return lin_int((0, 1), (n, 0), x)
        else:
            return lin_int((0, 1), (-n, 0), x)

    basis = np.vectorize(basis_1d, [np.double])


    stencil = SparseStencil()

    for p in positions:
        weight = basis(p, coarsening).prod()
        stencil.append(p, weight)

    return stencil


def ml_interpolation(fine_grid, coarse_grid):
    """Multilinear interpolation

    :rtype: Node
    """
    d = fine_grid.dimension()
    ones = np.ones(d)

    assert(fine_grid.dimension() == coarse_grid.dimension())

    return \
        operator.from_stencil(
            ml_interpolation_stencil(fine_grid, coarse_grid), fine_grid) * \
        operator.injection_interpolation(fine_grid, coarse_grid)

def fw_restriction_stencil(fine_grid, coarse_grid):
    """Return the stencil of a full weighting restriction.

    :rtype: SparseStencil
    """

    # compute the full weighting stencil from the interpolation
    interpolation = ml_interpolation_stencil(fine_grid, coarse_grid)

    s = 0
    for o, v in interpolation:
        s += v

    return (1/s) * interpolation.adjoint()


def fw_restriction(fine_grid, coarse_grid):
    """Return full weighting restriction operator.

    :rtype: Node
    """
    d = fine_grid.dimension()
    ones = np.ones(d)

    return \
        operator.injection_restriction(fine_grid, coarse_grid) * \
        operator.from_stencil(fw_restriction_stencil(fine_grid, coarse_grid),
        fine_grid)


def poisson_2d(grid, eps = 1.0):
    r"""The stencil of the discrete Poisson equation in 2D.

    This operator is the discrete version of the operator :math:`L` given by

    .. math::
       L u = -\left(
             \epsilon \frac{\partial^2 u}{\partial x^2} +
             \frac{\partial^2 u}{\partial y^2}
             \right)
             \,.

    Using finite differences leads to the discrete operator :math:`L_h`, which
    is computed by this function (see :ref:`poisson_equation`).

    :rtype: StencilNode
    """

    h0, h1 = grid.step_size()
    entries = [
        (( 0, -1), -1 / (h1*h1)),
        ((-1,  0), -1 / (h0*h0) * eps),
        (( 0,  0),  2 / (h0*h0) * eps + 2 / (h1*h1)),
        (( 1,  0), -1 / (h0*h0) * eps),
        (( 0,  1), -1 / (h1*h1))
    ]

    return operator.from_stencil(entries, grid)

def poisson_1d(grid):
    r"""The stencil of the discrete Poisson equation in 1D.

    :rtype: StencilNode
    """

    h, = grid.step_size()
    entries = [
        ((-1, ), -1 / (h*h)),
        (( 0, ),  2 / (h*h)),
        (( 1, ), -1 / (h*h)),
    ]

    return operator.from_stencil(entries, grid)


