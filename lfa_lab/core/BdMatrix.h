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

#ifndef LFA_BDMATRIX_H
#define LFA_BDMATRIX_H

#include "Common.h"

namespace lfa {

/** Block diagonal matrix storage. */
class BdMatrix {

  public:
    BdMatrix(int no_diag_blocks = 0, int block_rows = 0, int block_cols = 0);

    /** WARNING: Destroys all stored entries. */
    void resize(int no_diag_blocks = 0, int block_rows = 0, int block_cols = 0);
    void setZero();

    MatrixXcd& block(int i) { return m_diag_matrices[i]; }
    const MatrixXcd& block(int i) const { return m_diag_matrices[i]; }

    MatrixXcd full() const;

    int no_blocks() const { return m_diag_matrices.size(); }
    int rows() const { return m_block_rows * no_blocks(); }
    int cols() const { return m_block_cols * no_blocks(); }
    int block_rows() const { return m_block_rows; }
    int block_cols() const { return m_block_cols; }

    BdMatrix operator+ (const BdMatrix& rhs) const;
    BdMatrix operator* (const BdMatrix& rhs) const;
    friend BdMatrix operator* (complex<double> scalar, const BdMatrix& mat);

    bool dimensions_match(const BdMatrix& rhs) const;

    BdMatrix inverse() const;
    BdMatrix adjoint() const;

    double squaredNorm() const;
    double norm() const;

    double spectral_radius() const;
    double spectral_norm() const;

    VectorXcd eigenvalues() const;
  private:
    int m_block_rows;
    int m_block_cols;

    vector<MatrixXcd> m_diag_matrices;
};

ostream& operator<< (ostream& os, const BdMatrix& m);

}

#endif
