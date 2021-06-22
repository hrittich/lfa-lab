# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2021  Hannah Rittich
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

import numpy as np
from lfa_lab.dag import *
import lfa_lab.operator as operator

__all__ = 'smoothing_factor', 'h_ellipticity'

def smoothing_factor(op,
                     coarsening=None,
                     desired_resolution = None,
                     base_frequency = None):
    """Computes the smoothing factor.

    :param op: The operator to analyze.
    :type op: lfa_lab.dag.Node
    :param coarsening: The coarsening factor per dimension.
    :type coarsening: Tuple[int, ...]
    :param desired_resolution: The sampling resolution.
    :type desired_resolution: Tuple[int, ...]
    :param base_frequency: The lowest sampled frequency.
    :type base_frequency: Tuple[float, ...]
    """

    fine = op.output_grid
    if coarsening is None:
        coarsening = (2,) * fine.dimension()
    coarse = fine.coarse(coarsening)

    P = operator.hp_filter(fine, coarse)

    return (P * op).symbol().spectral_radius()

def h_ellipticity(op,
                  coarsening=None,
                  desired_resolution = None,
                  base_frequency = None):
    """Computes the :math:`h`-ellipticity.

    :param op: The operator to analyze.
    :type op: lfa_lab.dag.Node
    :param coarsening: The coarsening factor per dimension.
    :type coarsening: Tuple[int, ...]
    :param desired_resolution: The sampling resolution.
    :type desired_resolution: Tuple[int, ...]
    :param base_frequency: The lowest sampled frequency.
    :type base_frequency: Tuple[float, ...]
    """

    fine = op.output_grid
    if coarsening is None:
        coarsening = (2,) * fine.dimension()
    coarse = fine.coarse(coarsening)

    low_mode_count = np.prod(np.array(coarsening))

    if base_frequency is None:
        base_frequency = (0,) * fine.dimension()

    P = operator.hp_filter(fine, coarse)
    filtered_ews = (P * op).symbol(
        desired_resolution=desired_resolution,
        base_frequency=base_frequency).eigenvalues()

    # Since we are filtering out the low mode by projecting them to zero we
    # need to ignore the corresponding zero eigenvalues.
    assert (len(filtered_ews) % low_mode_count) == 0
    ignore_count = len(filtered_ews) // low_mode_count

    min_non_projected_ew = np.sort(np.abs(filtered_ews))[ignore_count]

    # Compute the denominator
    r = op.symbol(desired_resolution=desired_resolution,
                  base_frequency=base_frequency).spectral_radius()

    return min_non_projected_ew / r
       


