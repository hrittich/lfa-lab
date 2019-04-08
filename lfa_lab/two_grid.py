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

def coarse_grid_correction( \
        operator,
        coarse_operator,
        interpolation,
        restriction,
        coarse_error=None):
    """The error propagator of a coarse grid correction.

    See :ref:`error_coarse_grid_correction`.

    :param Node operator: The operator of the linear system.
    :param Node coarse_operator: The operator of the coarse linear system.
    :param Node interpolation: The interpolation operator.
    :param Node restriction: The restriction operator.
    :param Node coarse_error:
        The coarse error propagation operation. This entry is optional. It can
        be specified when the coarse grid system is solved inexactly.
    :return: The error propagator of the coarse grid correction.
    :rtype: Node
    """
    grid = operator.output_grid
    coarse_grid = coarse_operator.output_grid
    I  = operator.matching_identity()
    Ic = coarse_operator.matching_identity()
    P = interpolation
    R = restriction
    L = operator
    Lc = coarse_operator
    Ec = coarse_error

    if coarse_error is None:
        Ec = Ic.matching_zero()

    return (I - P * (Ic - Ec) * Lc.inverse() * R * L)

def galerkin_coarsening( \
        operator,
        interpolation,
        restriction):
    r"""The error propagator of the Galerkin coarse grid approximation (GCA).

    The GCA is defined as

    .. math:: L_c = R L P \,.

    :param Node operator: The fine grid operator :math:`L`.
    :param Node interpolation: The interpolation operator :math:`P`.
    :param Node restriction: The restriction operator :math:`R`.
    :return: The error propagation operator.
    :rtype: Node
    """

    L = operator
    P = interpolation
    R = restriction

    return R * L * P


def two_grid(pre_smoother,
             post_smoother,
             coarse_grid_correction):
    """The error propagator of a two-grid method."""
    S1 = pre_smoother
    S2 = post_smoother
    CGC = coarse_grid_correction

    return S2 * CGC * S1

