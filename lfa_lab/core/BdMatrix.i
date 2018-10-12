/* vim: set filetype=cpp: */
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

%feature("autodoc", "Block diagonal matrix.") BdMatrix;
%feature("autodoc", "The full matrix.") BdMatrix::full;
%feature("autodoc", "The i-th block of the matrix.") BdMatrix::block;
class BdMatrix {
  public:
    MatrixXcd full() const;

    MatrixXcd block(int i);

    int no_blocks() const;
    int rows() const;
    int cols() const;
    int block_rows() const;
    int block_cols() const;
};


