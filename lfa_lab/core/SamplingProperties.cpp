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

#include "SamplingProperties.h"
#include "MathUtil.h"

namespace lfa {

  SamplingProperties::SamplingProperties(ArrayFi finest_resolution, Grid grid)
    : m_finest_resolution(finest_resolution)
  {
    // compute the default base frequency
    m_base_frequency =
      (2.0 * pi / (grid.finestStepSize()
                   * m_finest_resolution.cast<double>())) / 2.0;
  }

  SamplingProperties::SamplingProperties(ArrayFi finest_resolution, ArrayFd base_frequency)
    : m_finest_resolution(finest_resolution),
    m_base_frequency(base_frequency)
  {

  }

}

