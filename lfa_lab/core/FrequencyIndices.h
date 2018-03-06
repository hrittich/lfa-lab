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

#ifndef LFA_FREQUENCY_INDICES_H
#define LFA_FREQUENCY_INDICES_H

#include "Common.h"
#include "NdRange.h"
#include "HarmonicIndices.h"

namespace lfa {

class HarmonicClusters;

/**
 * Space of all frequency indices.
 */
class FrequencyIndices
{
  public:
    /**
     * @param freq_per_dim Number of frequencies for every dimension.
     */
    FrequencyIndices(ArrayFi freq_per_dim = ArrayFi::Zero(0));

    // convert the harmonic coordinate to an index
    //int toLinear(ArrayFi coord) { return m_idx(coord); }

    /** Split into subsets given by T.
     * @param nsubsets Number of subsets per dimension
     * */
    vector<HarmonicIndices> split(ArrayFi nsubsets);

    vector<HarmonicIndices> splitIntoSizeOf(ArrayFi subset_freq_per_dim);

    int dimension() const { return m_freq_per_dim.rows(); }

    vector<ArrayFi> elements() const;

    vector<ArrayFd> frequencies(Frequency base, ArrayFd step_size) const;

    vector<ArrayFd> shifts(ArrayFd step_size) const;

    int indexOf(ArrayFi coord) {
      NdRange g(m_freq_per_dim); return g.indexOf(coord);
    }

    size_t size() { NdRange g(m_freq_per_dim); return g.size(); }
  private:
    // the number of frequencies per dimension
    ArrayFi m_freq_per_dim;
};

}

#endif
