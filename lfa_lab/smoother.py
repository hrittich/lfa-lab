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

from dag import *

from util import *
from operator import *

__all__ = [
    'jacobi',
    'gs_lex',
    'rb_jacobi'
]

def jacobi(op, weight=1.0):
    """The Jacobi smoother.

    Given by

    .. math:: I - \omega D^{-1} A \,.

    where :math:`D` is the diagonal part of :math:`A`.

    :param StencilNode op: The original operator :math:`A`.
    :param weight: The weight :math:`\omega`.
    """
    grid = op.grid

    A = op
    I = identity(grid)
    D = A.diag()

    return (I - weight * D.inverse() * A)

def gs_lex(op):
    """The Gauss-Seidel lexicographic smoother.

    Given by

    .. math:: I - (D + L)^{-1} A

    where :math:`D` is the diagonal part and
    :math:`L` the strictly lower triangular part of
    :math:`A`.

    :param StencilNode op: The original operator :math:`A`.
    """
    grid = op.grid

    A = op
    I = identity(grid)
    D = op.diag()
    L = op.lower()

    return (I - (D + L).inverse() * A)

def rb_jacobi(stencil, weight=1.0):
    r"""The red-black Jacobi method.

    :param StencilNode stencil: The original operator.
    :param double weight:
        The weight of the Jacobi methods, see
        :py:func:`lfa_lab.smoother.jacobi`.
    """
    grid = stencil.grid

    period = grid.dimension() * (2,)
    red_filter = NdArray(shape=period)
    black_filter = NdArray(shape=period)

    for i in NdRange(period):
        i = np.array(i)
        if is_even(i.sum()):
            red_filter[i] = IdentityNode(grid)
            black_filter[i] = ZeroNode(grid)
        else:
            red_filter[i] = ZeroNode(grid)
            black_filter[i] = IdentityNode(grid)

    red_filter = BlockNode(red_filter._entries)
    black_filter = BlockNode(black_filter._entries)

    J = jacobi(stencil, weight)
    return (red_filter + black_filter * J) * \
           (black_filter + red_filter * J)


