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

"""This module contains methods for construction block smoothers."""



from .util import *
from .core import *
import numpy as np
from .stencil import *
from .dag import *
from .operator import *

def _block_diag_stencil(stencil, block_size):

    indices = NdRange(block_size)
    stencils = NdArray(shape=block_size)

    for i in indices:
        i = np.array(i, np.int)
        s = stencil.filter( \
                lambda o, v: indices.inRange(i+o))
        stencils[i] = s

    return PeriodicStencil(stencils._entries)


def block_jacobi(op, block_size, weight = 1.0):
    """Returns the operator of the block jacobi method.

    :param StencilNode op: The original operator.
    :param tuple block_size: A tuple containing the block size per dimension.
    :param weight: The weight applied to the correction.
    """
    block_size = np.array(block_size)

    grid = op.grid
    diag_stencil = _block_diag_stencil(op.stencil, block_size)

    A = op
    I = identity(grid)
    D = from_periodic_stencil(diag_stencil, op.grid)

    return (I - weight * D.inverse() * A)

def rb_block_jacobi(op, block_size, weight = 1.0):
    """The red-black block Jacobi method.

    :param StencilNode op: The original operator.
    :param tuple block_size: A tuple containing the block size per dimension.
    :param weight: The weight applied to the correction.
    :return: The error propagation operator of the method.
    :rtype: Node
    """
    grid = stencil.grid

    block_size = np.array(block_size)

    red_filter = NdArray(shape=2*block_size)
    black_filter = NdArray(shape=2*block_size)

    for i in NdRange(2*block_size):
        i = np.array(i)
        p = i // block_size  # the number of the block
        if is_even(p.sum()):
            # red
            red_filter[i] = IdentityNode(grid)
            black_filter[i] = ZeroNode(grid)
        else:
            # black
            red_filter[i] = ZeroNode(grid)
            black_filter[i] = IdentityNode(grid)

    red_filter = BlockNode(red_filter._entries)
    black_filter = BlockNode(black_filter._entries)

    J = block_jacobi(stencil, block_size, weight)
    return (red_filter + black_filter * J) * \
           (black_filter + red_filter * J)

