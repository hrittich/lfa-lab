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

#include <gtest/gtest.h>

#include "BdMatrix.h"

using namespace lfa;

TEST(BdMatrix, Full)
{
    MatrixXcd A(2, 2);
    A << 1, 2,
         3, 4;

    MatrixXcd B(2,2);
    B << 5, 6,
         7, 8;

    MatrixXcd C(2,2);
    C <<  9, 10,
         11, 12;

    MatrixXcd D(6, 6);
    D << 1, 2, 0, 0, 0, 0,
         3, 4, 0, 0, 0, 0,
         0, 0, 5, 6, 0, 0,
         0, 0, 7, 8, 0, 0,
         0, 0, 0, 0, 9,10,
         0, 0, 0, 0,11,12;

    BdMatrix M(3, 2,2);
    M.set_block(0, A);
    M.set_block(1, B);
    M.set_block(2, C);

    EXPECT_LE((M.full() - D).norm(), 1e-10);
}

/** Lexicographic less for complex numbers. */
bool cmplx_lex_less(complex<double> a, complex<double> b)
{
    return (real(a) < real(b)) ||
        (real(a) == real(b) && imag(a) < imag(b));
}

TEST(BdMatrix, eigenvalues)
{
    MatrixXcd A(2, 2);
    A << 3, 0,
         0, 1;
    MatrixXcd B(2, 2);
    B << 4, 0,
         0, 2;

    MatrixXcd C(2, 2);
    C << 7, 0,
         0, 7;

    BdMatrix M(3, 2,2);
    M.set_block(0, A);
    M.set_block(1, B);
    M.set_block(2, C);

    VectorXcd ews = M.eigenvalues();


    VectorXcd expect(6);
    expect << 1, 2, 3, 4, 7, 7;

    // sort ews
    std::sort(ews.data(), ews.data()+ews.size(), cmplx_lex_less);
    std::sort(expect.data(), expect.data()+expect.size(), cmplx_lex_less);

    EXPECT_LE((ews - expect).norm(), 1e-10);
}


