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

#include "BlockStencil.h"

namespace lfa {

  BlockStencil BlockStencil::upper() const
  {
    BlockStencil U(shape());

    for (ConstIterator it(*this); it; ++it)
    {
      U(it.pos()) = it.value().upper();
    }

    return U;
  }


  BlockStencil BlockStencil::diag() const
  {
    BlockStencil D(shape());

    for (ConstIterator it(*this); it; ++it)
    {
      D(it.pos()) = it.value().diag();
    }

    return D;
  }

  BlockStencil BlockStencil::lower() const
  {
    BlockStencil L(shape());

    for (ConstIterator it(*this); it; ++it)
    {
      L(it.pos()) = it.value().lower();
    }

    return L;
  }

  BlockStencil BlockStencil::adjoint()
  {
    MultiArray< vector<StencilElement> >
      aux(ArrayFi::Zero(dimension()),
          shape() - ArrayFi::Ones(dimension()));

    for (Iterator block_it(*this); block_it; ++block_it)
    {
      ArrayFi a_offset = block_it.pos();
      DenseStencil a = block_it.value();

      // loop through the elements of the stencil a this position
      for (DenseStencil::ConstIterator it(a); it; ++it)
      {
        StencilElement adj_el;

        adj_el.offset = -it.pos(); // reverse the direction
        adj_el.value = it.value(); // ToDo Conjugate...

        // find the endpoint
        ArrayFi other_pos = mod(a_offset + it.pos(), shape());
        aux(other_pos).push_back(adj_el);
      }
    }

    // convert lists to block stencil
    BlockStencil B(shape());

    for (MultiArray< vector<StencilElement> >::Iterator it(aux); it; ++it)
    {
      B(it.pos()).setFromList(it.value());
    }

    return B;
  }

  BlockStencil BlockStencil::coarse()
  {
    // construct compatible size
    ArrayFi new_shape(dimension());
    for (int i = 0; i < dimension(); ++i) {
      assert( shape()(i) % 2 == 0 || shape()(i) == 1);

      if (shape()(i) > 1) {
        new_shape(i) = shape()(i) / 2;
      } else {
        new_shape(i) = shape()(i);
      }
    }

    BlockStencil B(new_shape);

    for (BlockStencil::Iterator it(B); it; ++it) {
      B(it.pos()) = (*this)(2 * it.pos()).coarse( 2*ArrayFi::Ones(dimension()) );
    }

    return B;
  }

  DenseStencil BlockStencil::multiplyRight(const DenseStencil& a,
                                           ArrayFi a_offset) const
  {
    DenseStencil c;

    bool min_max_uninit = true;

    ArrayFi start;
    ArrayFi end;

    // compute the size of the current stencil entry
    for (DenseStencil::ConstIterator it(a); it; ++it) {

      // the stencil at (block_it.pos() + it.pos())
      ArrayFi other_pos = mod(a_offset + it.pos(), shape());
      const DenseStencil& b = (*this)(other_pos);

      ArrayFi cur_start = it.pos() + b.startIndex();
      ArrayFi cur_end = it.pos() + b.endIndex();

      if (min_max_uninit) {
        start = cur_start;
        end = cur_end;

        min_max_uninit = false;
      } else {
        start = cur_start.min(start);
        end = cur_end.max(end);
      }
    }

    // resize and fill with zero
    c.resize(start, end);

    // compute the stencil coefficients
    for (DenseStencil::ConstIterator it(a); it; ++it) {

      // the stencil at (block_it.pos() + it.pos())
      ArrayFi other_pos = mod(a_offset + it.pos(), shape());
      const DenseStencil& b = (*this)(other_pos);

      double coeff_a = it.value();

      // compute coefficients
      for (DenseStencil::ConstIterator inner_it(b); inner_it; ++inner_it)
      {
        ArrayFi dest_pos = it.pos() + inner_it.pos();
        c(dest_pos) += coeff_a * inner_it.value();
      }
    }

    return c;
  }

  BlockStencil operator* (const BlockStencil& A, const BlockStencil& B)
  {
    BlockStencil C(A.shape()); // result

    // assume the block sizes are the same, so the block size will stay the
    // same
    assert( (A.shape() == B.shape()).all() );


    for (BlockStencil::Iterator block_it(C); block_it; ++block_it)
    {
      ArrayFi a_offset = block_it.pos();
      const DenseStencil& a = A(a_offset);

      C(block_it.pos()) = B.multiplyRight(a, a_offset);
    }

    return C;
  }

  BlockStencil operator* (double s, const BlockStencil& A)
  {
    BlockStencil B(A.shape());

    for (BlockStencil::ConstIterator it(A); it; ++it) {
      B(it.pos()) = s * it.value();
    }

    return B;
  }

  std::ostream& operator<< (std::ostream& os, const BlockStencil& A)
  {
    for (BlockStencil::ConstIterator it(A); it; ++it)
    {
      os << "=== " << it.pos().transpose() << " ===" << std::endl;
      os << it.value() << std::endl;
    }

    return os;
  }

}
