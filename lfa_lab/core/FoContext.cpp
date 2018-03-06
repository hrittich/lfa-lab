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

#include "FoContext.h"
#include "MathUtil.h"

namespace lfa {

  FoContext::FoContext(int dimension, ArrayFd step_size)
    : m_dimension(dimension),
    m_step_size(step_size)
  {
    if (m_step_size.rows() == 0) {
      m_step_size = 1.0 / 64 * ArrayFd::Ones(dimension);
    }

    if (m_step_size.rows() != dimension) {
      throw logic_error("Step size has invalid number of rows.");
    }
  }

}

