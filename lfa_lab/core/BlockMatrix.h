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

#ifndef LFA_BLOCK_MATRIX_H
#define LFA_BLOCK_MATRIX_H

#include <Common.h>

namespace lfa {

/** A matrix of matrices. */
class BlockMatrix {
    public:

        /**
         * @param block_rows Number of rows per block.
         * @param nrows Number of block rows.
         */
        BlockMatrix(int block_rows = 0, int block_cols = 0, int nrows = 0, int ncols = 0);

        BlockMatrix(int block_rows, int block_cols, int nrows, int ncols, const MatrixXcd& store);

        Block<MatrixXcd> operator() (int row, int col) {
            return m_store.block(row * m_block_rows,
                                 col * m_block_cols,
                                 m_block_rows,
                                 m_block_cols);
        }

        friend BlockMatrix operator* (complex<double> s, const BlockMatrix& A);
        friend BlockMatrix operator* (const BlockMatrix& A, const BlockMatrix& B);

        BlockMatrix operator+ (const BlockMatrix& B) {
            assert(m_nrows == B.m_nrows);
            assert(m_ncols == B.m_ncols);
            assert(m_block_rows == B.m_block_rows);
            assert(m_block_cols == B.m_block_cols);

            return BlockMatrix(m_block_rows, m_block_cols,
                               m_nrows, m_ncols,
                               m_store + B.m_store);
        }

        const MatrixXcd& toMatrix() { return m_store; }

        BlockMatrix inverse() const;

    private:
        int m_block_rows;
        int m_block_cols;
        int m_nrows;
        int m_ncols;

        MatrixXcd m_store;
};

inline BlockMatrix operator* (complex<double> s, const BlockMatrix& A) {
    return BlockMatrix(A.m_block_rows, A.m_block_cols,
                       A.m_nrows, A.m_ncols,
                       s * A.m_store);
}

inline BlockMatrix operator* (const BlockMatrix& A, const BlockMatrix& B)
{
    assert(A.m_ncols == B.m_nrows);
    assert(A.m_block_cols == B.m_block_rows);

    return BlockMatrix(A.m_block_rows, B.m_block_cols,
                       A.m_nrows, B.m_ncols,
                       A.m_store*B.m_store);
}

}

#endif
