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

#ifndef LFA_HARMONIC_INDICES_H
#define LFA_HARMONIC_INDICES_H

#include "Common.h"

namespace lfa {

  class HarmonicIndices
  {
    public:
      HarmonicIndices(ArrayFi freq_per_dim = ArrayFi::Zero(0),
                      ArrayFi multiplier = ArrayFi::Zero(0),
                      ArrayFi shift = ArrayFi::Zero(0));

      /** Return the elements of the set as a vector. */
      vector<ArrayFi> elements() const;

      /** Return the indices of the elements in this set. (With respect to
       * the super set.)
       */
      vector<int> indices() const;

      /** Give the base frequency for this subset when given the base
       * frequency for the superset. */
      ArrayFd baseFrequency(ArrayFd base, ArrayFd step_size) const;

      bool operator== (const HarmonicIndices& other) const;
      bool operator!= (const HarmonicIndices& other) const {
        return ! ((*this) == other);
      }
    private:
      // number of frequencies in this set
      ArrayFi m_freq_per_dim;

      // number of subsets in the original set
      ArrayFi m_multiplier;

      // which subset is this?
      ArrayFi m_shift;
  };

}

#endif
