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

#ifndef LFA_HARMONIC_CLUSTERS_H
#define LFA_HARMONIC_CLUSTERS_H

#include "Common.h"
#include "NdRange.h"
#include "FrequencyIndices.h"
#include "FoContext.h"

namespace lfa {

  /** Indexing for the harmonics. */
  class HarmonicClusters {
    public:
      explicit HarmonicClusters(ArrayFi base_shape = ArrayFi::Zero(0),
                                ArrayFi cluster_shape = ArrayFi::Zero(0));

      bool operator== (const HarmonicClusters& other) const;
      bool operator!= (const HarmonicClusters& other) const {
        return !((*this) == other);
      }

      /** Shape of the global(!) index space.
       * In most cases this is NOT what you want!
       */
      ArrayFi shape() const;

      /** Merge clusters into larger ones. */
      HarmonicClusters mergeCluster(ArrayFi factor) const;

      NdRange baseIndices() const { return NdRange(m_base_shape); }
      NdRange clusterIndices() const { return NdRange(m_cluster_shape); }
      NdRange globalIndices() const { return NdRange(shape()); }

      int clusterSize() const;
      int baseSize() const;
      int size() const;

      /** Two harmonic clusters are compatible if they have the same base
       * shape. */
      bool isCompatibleTo(const HarmonicClusters& other) const;

      int dimension() const { return m_base_shape.rows(); }

      ArrayFi globalIndex(ArrayFi base_index, ArrayFi cluster_index) const;

      ArrayFi baseIndex(ArrayFi global_index) const;
      ArrayFi clusterIndex(ArrayFi global_index) const;

      /** Convert a coordinate in this coordinate system to a coordinate in
       * the result_clusters coordinate system. */
      void convert(ArrayFi& result_base_index,
                   ArrayFi& result_cluster_index,
                   HarmonicClusters result_clusters,
                   ArrayFi base_index,
                   ArrayFi cluster_index) const;

      /** HarmonicClusters with smallest possible cluster shape that is a
       * multiple of this and other. */
      HarmonicClusters minContainer(const HarmonicClusters& other) const;

      ArrayFi expansionFactor(HarmonicClusters& expanded) const;

      ArrayFi clusterShape() const { return m_cluster_shape; }
    private:
      ArrayFi m_base_shape;
      ArrayFi m_cluster_shape;
  };

}

#endif
