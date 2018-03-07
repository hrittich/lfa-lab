/*
  vim: set filetype=cpp:

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


template <typename T>
class MatrixContainer
{
  public:
    MatrixContainer(int nrows = 0, int ncols = 0);

    T& operator() (int row, int col);

    void resize(int rows, int cols);

    int rows() const;
    int cols() const;

    bool empty() const;

    %extend {
      T __getitem__(int i, int j) {
        return (*$self)(i, j);
      }

      void __setitem__(ArrayFi p, T value) {
        if (p.rows() != 2) {
          throw std::logic_error("Invalid index length.");
        }
        (*$self)(p[0], p[1]) = value;
      }
    }
};


