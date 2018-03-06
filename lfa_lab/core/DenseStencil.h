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

#ifndef LFA_DENSE_STENCIL_H
#define LFA_DENSE_STENCIL_H

#include "Common.h"
#include "MultiArray.h"
#include "StencilElement.h"

namespace lfa {

  /** Storage for a dense n-dimensional stencil. */
  class DenseStencil : public MultiArray<double>
  {
    public:
      /** Create a stencil. The dimension is given by the length of the two
       * vectors.
       * The stencil is initially filled with zeros. */
      DenseStencil(ArrayFi start = ArrayFi::Zero(1),
                   ArrayFi end = ArrayFi::Zero(1));

      /** Sets the stencil to the given size and fill it with zeros. */
      void resize(ArrayFi start, ArrayFi end);

      void setZero(int dim);
      void setIdentity(int dim);


      typedef vector<StencilElement> ElementList;
      /** Sets the stencil from a list of elements. The stencil will be
       * resized appropriately. */
      void setFromList(const ElementList& list);

      /** A stencil representing the diagonal element. */
      DenseStencil diag() const;

      /** A stencil representing the strictly lower diagonal elements. */
      DenseStencil lower() const;

      /** A stencil representing the strictly upper diagonal elements. */
      DenseStencil upper() const;

      DenseStencil coarse(const ArrayFi& space);
  };

  /** Multiply two stencils. */
  DenseStencil operator* (const DenseStencil& s, const DenseStencil& t);

  /** Multiply a DenseStencil by a scalar value  */
  DenseStencil operator* (double s, const DenseStencil& t);

  /** Computes the galerkin coarse grid stencil. */
  DenseStencil galerkin_stencil(
          const DenseStencil& R,
          const DenseStencil& L,
          const DenseStencil& P);

}

#endif
