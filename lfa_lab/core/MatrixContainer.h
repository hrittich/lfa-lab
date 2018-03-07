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

#ifndef LFA_MATRIX_CONTAINER_H
#define LFA_MATRIX_CONTAINER_H

#include "Common.h"

namespace lfa {

  template <typename T>
  class MatrixContainer
  {
    public:
      MatrixContainer(int nrows = 0, int ncols = 0)
        :  m_store(nrows, vector<T>(ncols)),
        m_nrows(nrows),
        m_ncols(ncols)
      {}

      template <typename R>
      MatrixContainer(const MatrixContainer<R>& rhs)
      :  m_nrows(0), m_ncols(0)
      {
        resize(rhs.rows(), rhs.cols());

        for (int i = 0; i < rows(); ++i) {
          for (int j = 0; j < cols(); ++j) {
            (*this)(i,j) = rhs(i,j);
          }
        }
      }

      T& operator() (int row, int col) {
        return m_store[row][col];
      }

      const T& operator() (int row, int col) const {
        return m_store[row][col];
      }

      void resize(int rows, int cols)
      {
        m_nrows = rows;
        m_ncols = cols;

        m_store.resize(rows);
        for (int i = 0; i < rows; ++i)
        {
          m_store[i].resize(cols);
        }
      }

      int rows() const { return m_nrows; }
      int cols() const { return m_ncols; }

      bool empty() const {
        return (rows() == 0 || cols() == 0);
      }
    private:
      vector< vector<T> > m_store;

      int m_nrows;
      int m_ncols;
  };

}

#endif
