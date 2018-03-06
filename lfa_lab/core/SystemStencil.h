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

#ifndef LFA_SYSTEM_STENCIL_H
#define LFA_SYSTEM_STENCIL_H

#include "Common.h"
#include "DenseStencil.h"
#include <algorithm>
#include "MatrixContainer.h"

namespace lfa {

  /** Storage for a rows x cols system of stencil equation. */
  class SystemStencil : public MatrixContainer<DenseStencil>
  {
    public:
      /** Initialize a scalar stencil system. */
      SystemStencil(const DenseStencil& s);
      /** Create an equation system consisting of the given number of cols
       * and rows. */
      SystemStencil(int rows = 0, int cols = 0);

      /** Returns the maximum of cols and rows. */
      int grids() const { return std::max(rows(), cols()); }

      /** The dimension of the stencils. All stencils should have the same
       * dimension. */
      int dimension() const {
        assert(!empty());
        return (*this)(0,0).dimension();
      }
  };

  SystemStencil galerkin_stencil(
      const SystemStencil& R,
      const SystemStencil& A,
      const SystemStencil& P);

}

#endif
