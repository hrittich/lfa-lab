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

#ifndef LFA_DISCRETE_DOMAIN_H
#define LFA_DISCRETE_DOMAIN_H

#include "Common.h"

#include "SplitFrequencyDomain.h"
#include "SamplingProperties.h"

namespace lfa {

  /** Discrete subset of a SplitFrequencyDomain.
   * This is used for sampling frequency operators.
   */
  class DiscreteDomain
  {
    public:
      DiscreteDomain(
          const SplitFrequencyDomain& domain,
          const SamplingProperties& conf);

      HarmonicClusters harmonics() const;

      ArrayFd frequency(ArrayFi base_index, ArrayFi cluster_index) const;
      ArrayFd frequency(ArrayFi global_index) const;

      int dimension() const { return m_domain.dimension(); }

      /** The resolution on the current grid. */
      ArrayFi resolution() const;
      /** Step size of the current grid. */
      ArrayFd step_size() const;

      Grid grid() const { return m_domain.grid(); }
    private:
      SplitFrequencyDomain m_domain;
      SamplingProperties m_conf;

  };

}

#endif
