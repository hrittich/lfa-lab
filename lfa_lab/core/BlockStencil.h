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

#ifndef LFA_BLOCK_STENCIL_H
#define LFA_BLOCK_STENCIL_H

#include "Common.h"
#include "DenseStencil.h"
#include "MultiArray.h"
#include "MathUtil.h"

namespace lfa {

  class BlockStencil : public MultiArray<DenseStencil>
  {
    public:
      BlockStencil(const BlockStencil& rhs)
        : MultiArray<DenseStencil>(rhs)
      { }

      BlockStencil(ArrayFi shape = ArrayFi::Ones(1))
        : MultiArray(
            ArrayFi::Zero(shape.rows()),
            shape - ArrayFi::Ones(shape.rows()))
      {
      }

      void reshape(ArrayFi shape) {
        resize(ArrayFi::Zero(shape.rows()),
            shape - ArrayFi::Ones(shape.rows()));
      }

      ArrayFi shape() const {
        return ArrayFi(endIndex() - startIndex() +
            ArrayFi::Ones(startIndex().rows()));
      }

      /** Same number of entries per dimension. */
      inline bool isSquare() const {
        ArrayFi sh = shape();
        for (int d = 0; d < dimension()-1; ++d) {
          if (sh(d) != sh(d+1))
            return false;
        }
        return true;
      }

      inline int frequencyRange() const {
        assert(isSquare());
        assert(isPow2( shape()(0) ));

        return log2_rn( shape()(0) );
      }

      BlockStencil upper() const;

      /** The substencil representing the diagonal elements only. */
      BlockStencil diag() const;

      /** The substencil representing the sub diagonal elements. */
      BlockStencil lower() const;

      /** Return the stencil for the adjoint of this operator. */
      BlockStencil adjoint();

      /** Compute a coarse representation .*/
      BlockStencil coarse();

      /** Multiply this blockstencil to the right hand side of a stencil.
       *
       * @param a The stencil on the left.
       * @param a_offset The position where the stencil a is positioned in
       *                 the block stencil.
       */
      DenseStencil multiplyRight(const DenseStencil& a, ArrayFi a_offset) const;

  };

  BlockStencil operator* (const BlockStencil& A, const BlockStencil& B);
  std::ostream& operator<< (std::ostream& os, const BlockStencil& A);

  BlockStencil operator* (double s, const BlockStencil& A);

}

#endif
