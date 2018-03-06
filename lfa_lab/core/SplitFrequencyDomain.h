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

#ifndef LFA_FREQUENCY_COUPLING_H
#define LFA_FREQUENCY_COUPLING_H

#include "Common.h"
#include "Grid.h"
#include "HarmonicClusters.h"

namespace lfa {

  class SplitFrequencyDomain {
    public:
      /** The domain of split frequencies.
       * If cluster_shape is omited then a cluster_shape of one is assumed.
       */
      explicit SplitFrequencyDomain(Grid grid = Grid(),
          ArrayFi cluster_shape = ArrayFi::Zero(0));

      const Grid& grid() const { return m_grid; }
      const ArrayFi& clusterShape() const { return m_cluster_shape; }

      /** Returns the harmonic cluster.
       *
       * @param resolution The number of sampling points on the current grid.
       */
      HarmonicClusters harmonics(ArrayFi resolution) const;

      /** Computes the resolution on the current grid from the resolution on
       * the finest grid.
       */
      ArrayFi resolution(ArrayFi finest_resolution) const;

      bool operator== (const SplitFrequencyDomain& other) const {
        return (m_grid == other.m_grid
            && (m_cluster_shape == other.m_cluster_shape).all());
      }
      bool operator!= (const SplitFrequencyDomain& other) const {
        return !((*this) == other);
      }

      bool isValid();

      ArrayFd step_size() const { return m_grid.step_size(); }

      int dimension() const { return m_grid.dimension(); }


      /** Returns if this and other can be the domains of the same
       * block Toeplitz operator. */
      bool isCompatibleTo(const SplitFrequencyDomain& other);

      SplitFrequencyDomain expand(ArrayFi factor) const;

      /** Compute the least common coupling of this and other. */
      ArrayFi lcc(const SplitFrequencyDomain& other);
    private:
      Grid m_grid;
      ArrayFi m_cluster_shape;
  };

  ostream& operator<< (ostream& os, const SplitFrequencyDomain& domain);

}

#endif
