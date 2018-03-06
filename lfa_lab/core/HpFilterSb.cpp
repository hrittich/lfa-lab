/*
  LFA Lab - Library to simplify local Fourier analysis.
  Copyright (C) 2018  Hannah Rittich

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "HpFilterSb.h"

#include "MathUtil.h"
#include "DiscreteDomain.h"


namespace lfa {

  HpFilterSb::HpFilterSb(Grid fine_grid, Grid coarse_grid)
    :  m_grid(fine_grid),
    m_coarsing_factor(fine_grid.coarsening_factor(coarse_grid))
  { }

  FoProperties HpFilterSb::properties()
  {
    SplitFrequencyDomain coupling(m_grid, ArrayFi::Ones(m_grid.dimension()));
    return FoProperties(coupling, coupling);
  }

  Symbol HpFilterSb::generate(const SamplingProperties& conf)
  {
    /*
       x = low frequency

       --------
       |x    x|
       |      |
       |      |
       |x    x|
       --------
       0       2\pi/h
       */

    int d = m_grid.dimension();

    // the size of the low frequency block
    ArrayFd block_size = pi / (m_grid.step_size() * m_coarsing_factor.cast<double>());

    ArrayFd lower_bound = block_size;
    ArrayFd upper_bound = (2 * m_coarsing_factor.cast<double>()
        - ArrayFd::Ones(d)) * block_size;


    ArrayFi zero = ArrayFi::Zero(d);

    DiscreteDomain domain(properties().input(), conf);
    HarmonicClusters cluster = domain.harmonics();
    Symbol result(cluster, cluster);

    NdRange bases = cluster.baseIndices();
    for (NdRange::iterator b = bases.begin(); b != bases.end(); ++b)
    {
      ArrayFd freq = domain.frequency(*b, zero);

      bool is_high = (freq > lower_bound && freq <= upper_bound).any();

      result.ref(*b, zero, zero) = (is_high ? 1 : 0);
    }

    return result;
  }


}
