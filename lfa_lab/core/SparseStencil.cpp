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

#include "SparseStencil.h"

namespace lfa {

SparseStencil::SparseStencil()
{ }

SparseStencil::SparseStencil(const DenseStencil& other)
{
    for (DenseStencil::ConstIterator p(other); p; ++p) {
        this->append(p.pos(), *p);
    }
}

int SparseStencil::dimension() const {
    if (m_elements.empty()) {
        return 0;
    }
    return m_elements.front().offset.rows();
}

void SparseStencil::append(ArrayFi offset, complex<double> value)
{
    if (dimension() != 0 && offset.rows() != dimension()) {
        std::stringstream ss;
        ss << "Inconsistent dimensions. Stencil has dimension "
           << dimension() << " but offset has dimension "
           << offset.rows() << ".";
        throw logic_error(ss.str());
    }
    m_elements.push_back(StencilElement(offset, value));
}



}

