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

#ifndef LFA_N_MATRIX_VIEW
#define LFA_N_MATRIX_VIEW

#include "Common.h"
#include "VectorizedIndex.h"

namespace lfa {

  /** Matrix view with n-dimensional indexing. */
  template <typename T>
    class NMatrixView
    {
      public:
        NMatrixView(T& matrix, ArrayFi shape, bool initToZero = false)
          : matrix(matrix),
            idx(ArrayFi::Zero(shape.rows()),
                shape - ArrayFi::Ones(shape.rows()))
      {
        if (initToZero) {
          matrix = T::Zero(idx.elements(), idx.elements());
        } else {
          assert(idx.elements() == matrix.rows());
          assert(idx.elements() == matrix.cols());
        }
      }

        typename T::Scalar& operator() (const ArrayFi& pos_i, const ArrayFi& pos_j)
        {
          return matrix(idx(pos_i), idx(pos_j));
        }

        typename T::Scalar operator() (const ArrayFi& pos_i, const ArrayFi& pos_j) const
        {
          return matrix(idx(pos_i), idx(pos_j));
        }

        const VectorizedIndex& indices() { return idx; }

      private:
        T& matrix;
        VectorizedIndex idx;
    };

}

#endif
