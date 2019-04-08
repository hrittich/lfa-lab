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

from .dag import *
from .stencil import *
from .util import *

def identity(grid):
    """The identity operator.

    :param Grid grid: The grid where the operator is defined.
    """
    return IdentityNode(grid)

def zero(grid):
    """The zero operator.

    :param Grid grid: The grid where the operator is defined.
    """
    return ZeroNode(grid)

def from_stencil(stencil, grid):
    """Create an operator from a stencil.

    See also :ref:`defining_stencil_operators`.

    :param stencil: The stencil that should be transformed into an operator.
      This parameter has to be assigned to either a SparseStencil or a list of
      tuples, where each tuple consists of an offset and a value.
    :param Grid grid: The grid where the operator is defined.
    """

    if not isinstance(stencil, SparseStencil):
        stencil = SparseStencil(entries = stencil)

    return StencilNode(stencil, grid)

def hp_filter(fine_grid, coarse_grid):
    """Create a high-pass filter.

    :param Grid fine_grid: The grid corresponding to all frequencies.
    :param Grid coarse_grid: The grid of the low modes.
    """
    return HpFilterNode(fine_grid, coarse_grid)

def lp_filter(fine_grid, coarse_grid):
    """Create a low-pass filter.

    :param Grid fine_grid: The grid corresponding to all frequencies.
    :param Grid coarse_grid: The grid of the low modes.
    """
    I = IdentityNode(fine_grid)
    HP = HpFilterNode(fine_grid, coarse_grid)
    return I - HP

def injection_interpolation(fine_grid, coarse_grid):
    """Create an injection interpolation operator.

    :param Grid fine_grid: The codomain of the operator.
    :param Grid coarse_grid: The domain of the operator.
    """
    return FlatInterpolationNode(fine_grid, coarse_grid)

def injection_restriction(fine_grid, coarse_grid):
    """Create an injection restriction operator.

    :param Grid fine_grid: The domain of the operator.
    :param Grid coarse_grid: The codomain of the operator.
    """
    return FlatRestrictionNode(coarse_grid, fine_grid)

def from_periodic_stencil(stencils, grid):
    """Create an operator from a periodic stencil.

    :param NdArray stencils: An :math:`n`-D array of stencils.
    """
    def wrap_stencil(s):
        if not isinstance(s, SparseStencil):
            return SparseStencil(entries = s)
        else:
            return s

    if not isinstance(stencils, PeriodicStencil):
        assert(isinstance(stencils, NdArray))
        stencils = PeriodicStencil(stencils.map(wrap_stencil))

    return PeriodicStencilNode(stencils, grid)

def system(entries):
    """Given a list of lists, returns a system of operators."""

    return SystemNode(entries)


