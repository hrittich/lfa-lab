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

#include "ClusterSymbol.h"

namespace lfa {

ClusterSymbol::ClusterSymbol(ArrayFi output_shape, ArrayFi input_shape)
  : m_output_indices(output_shape),
    m_input_indices(input_shape),
    m_store(m_output_indices.size(),
            m_input_indices.size())
{

}

complex<double>& ClusterSymbol::operator() (ArrayFi row, ArrayFi col)
{
    return m_store(
        m_output_indices.indexOf(row),
        m_input_indices.indexOf(col));
}

void ClusterSymbol::setMatrix(MatrixXcd m)
{
    assert(m.rows() == m_output_indices.size());
    assert(m.cols() == m_input_indices.size());

    m_store = m;
}


}

