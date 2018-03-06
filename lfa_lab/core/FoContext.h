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

#ifndef LFA_FO_CONTEXT_H
#define LFA_FO_CONTEXT_H

#include "Common.h"

namespace lfa {

  /** The context of the frequency analysis, i.e., resolution, step-size, and
   * sampling base frequency. */
  class FoContext {
    public:
      FoContext(int dimension,
                ArrayFd step_size = ArrayFd::Zero(0));

      const ArrayFd& step_size() const { return m_step_size; }

      int dimension() const { return m_dimension; }

    private:
      int m_dimension;

      /** Step size of the finest grid. */
      ArrayFd m_step_size;
  };

}

#endif
