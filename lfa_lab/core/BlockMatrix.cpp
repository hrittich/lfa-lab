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

#include "BlockMatrix.h"
#include <Eigen/Dense>

namespace lfa {

BlockMatrix::BlockMatrix(int block_rows, int block_cols, int nrows, int ncols)
 :  m_block_rows(block_rows),
    m_block_cols(block_cols),
    m_nrows(nrows),
    m_ncols(ncols)
{
    m_store = MatrixXcd::Zero(block_rows * nrows, block_cols * ncols);
}

BlockMatrix::BlockMatrix(int block_rows, int block_cols, int nrows, int ncols, const MatrixXcd& store)
 :  m_block_rows(block_rows),
    m_block_cols(block_cols),
    m_nrows(nrows),
    m_ncols(ncols),
    m_store(store)
{

}

BlockMatrix BlockMatrix::inverse() const {
    return BlockMatrix(m_block_rows, m_block_cols,
                        m_nrows, m_ncols,
                        m_store.inverse());
}



}
