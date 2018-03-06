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

#ifndef SAMPLING_PROPERTIES_H
#define SAMPLING_PROPERTIES_H

#include "Common.h"
#include "Grid.h"

namespace lfa {

  class SamplingProperties {
    public:
      /** Set finest_resolution and compute the default base_frequency.
       * @param grid Is an arbitrary grid from the grid hierarchy.
       */
      SamplingProperties(ArrayFi finest_resolution, Grid grid);
      SamplingProperties(ArrayFi finest_resolution, ArrayFd base_frequency);

      /** Resolution on the finest grid. */
      const ArrayFi& finest_resolution() const { return m_finest_resolution; }
      const ArrayFd& base_frequency() const { return m_base_frequency; }
    private:
      ArrayFi m_finest_resolution; /// < The resolution on the finest grid.
      ArrayFd m_base_frequency;
  };

}

#endif // SAMPLING_PROPERTIES_H
