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

#include "DenseStencil.h"

#include "MathUtil.h"

namespace lfa {

DenseStencil::DenseStencil(ArrayFi start, ArrayFi end)
    : MultiArray(start, end)
{
    fill(0);
}

void DenseStencil::resize(ArrayFi start, ArrayFi end)
{
    MultiArray::resize(start, end);
    fill(0);
}

void DenseStencil::setZero(int dim)
{
    // resize implicitly sets the stencil to zero
    resize(ArrayFi::Zero(dim), ArrayFi::Zero(dim));
}

void DenseStencil::setIdentity(int dim)
{
    setZero(dim);
    operator() ( ArrayFi::Zero(dim) ) = 1.0;
}

void DenseStencil::setFromList(const ElementList& list)
{
    // compute min and max element

    ArrayFi start = list.front().offset;
    ArrayFi end = list.front().offset;

    for (ElementList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
        start = start.min(it->offset); // elementwise min/max
        end = end.max(it->offset);
    }

    resize(start, end);

    // sum the elements
    for (ElementList::const_iterator it = list.begin(); it != list.end(); ++it)
    {
        (*this)(it->offset) += it->value;
    }
}

DenseStencil DenseStencil::diag() const
{
    DenseStencil D(ArrayFi::Zero(dimension()),
              ArrayFi::Zero(dimension()));

    D(ArrayFi::Zero(dimension()))
        = (*this)(ArrayFi::Zero(dimension()));

    return D;
}

DenseStencil DenseStencil::lower() const
{
    ArrayFi new_end = endIndex();
    new_end(dimension()-1) = 0;

    DenseStencil L(startIndex(), new_end);

    for (DenseStencil::ConstIterator it(*this); it; ++it) {
        if (it.beforeCenter()) {
            L(it.pos()) = it.value();
        }
    }

    return L;
}

DenseStencil DenseStencil::upper() const
{
    ArrayFi new_start = startIndex();
    new_start(dimension()-1) = 0;

    DenseStencil L(new_start, endIndex());

    for (DenseStencil::ConstIterator it(*this); it; ++it) {
        if (it.behindCenter()) {
            L(it.pos()) = it.value();
        }
    }

    return L;
}



DenseStencil DenseStencil::coarse(const ArrayFi& space)
{
    ArrayFi new_start(space.rows());
    ArrayFi new_end(space.rows());
    for (int i = 0; i < new_start.rows(); ++i)
    {
        new_start(i) = div_rz( startIndex()(i), space(i) );
        new_end(i) = div_rz( endIndex()(i), space(i) );
    }

    DenseStencil r(new_start, new_end);

    for (DenseStencil::Iterator it(r); it; ++it)
    {
        it.value() = (*this)( space * it.pos().array() );
    }

    return r;
}

DenseStencil operator* (const DenseStencil& s, const DenseStencil& t)
{
    DenseStencil r(s.startIndex() + t.startIndex(), s.endIndex() + t.endIndex());

    for (DenseStencil::ConstIterator it_s(s); it_s; ++it_s) {
        for (DenseStencil::ConstIterator it_t(t); it_t; ++it_t) {
            r( it_s.pos() + it_t.pos() ) += it_s.value() * it_t.value();
        }
    }

    return r;
}

DenseStencil operator* (double s, const DenseStencil& t)
{
    DenseStencil r(t.startIndex(), t.endIndex());

    for (DenseStencil::ConstIterator it(t); it; ++it) {
        r(it.pos()) = s * it.value();
    }

    return r;
}


DenseStencil galerkin_stencil(
        const DenseStencil& R,
        const DenseStencil& L,
        const DenseStencil& P)
{
    return (R * L * P).coarse( ArrayFi::Ones(L.dimension()) * 2 );
}




}
