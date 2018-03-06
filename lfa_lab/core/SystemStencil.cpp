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

#include "SystemStencil.h"

namespace lfa {

SystemStencil::SystemStencil(const DenseStencil& s)
 :  MatrixContainer(1, 1)
{
    (*this)(0,0) = s;
}

SystemStencil::SystemStencil(int rows, int cols)
 :  MatrixContainer(rows, cols)
{
}


SystemStencil galerkin_stencil(
        const SystemStencil& R,
        const SystemStencil& A,
        const SystemStencil& P)
{
    SystemStencil Ac(A.rows(), A.cols());

    for (int i = 0; i < A.rows(); ++i)
    {
        for (int j = 0; j < A.cols(); ++j)
        {
            // Ac_{ij} = R_i * A_{ij} * P_j

            Ac(i,j) = galerkin_stencil(R(0,i), A(i,j), P(j,0));
        }
    }

    return Ac;
}


}
