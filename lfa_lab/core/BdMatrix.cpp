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

#include "BdMatrix.h"
#include "EigenSolver.h"
#include "ExEigenSolver.h"

#include <stdexcept>
#include <Eigen/Dense>

namespace lfa {

BdMatrix::BdMatrix(int no_diag_blocks, int block_rows, int block_cols)
  : m_block_rows(block_rows),
    m_block_cols(block_cols),
    m_diag_matrices(no_diag_blocks)
{
    assert(no_diag_blocks >= 0);
}

void BdMatrix::resize(int no_diag_blocks, int block_rows, int block_cols)
{
    m_diag_matrices.resize(no_diag_blocks);
    fill(m_diag_matrices.begin(), m_diag_matrices.end(), MatrixXcd());

    m_block_rows = block_rows;
    m_block_cols = block_cols;
}

void BdMatrix::setZero()
{
    for (int i = 0; i < no_blocks(); ++i) {
        block(i) = MatrixXcd::Zero(m_block_rows, m_block_cols);
    }
}

MatrixXcd BdMatrix::full() const
{
    MatrixXcd M = MatrixXcd::Zero(rows(), cols());

    for (int i = 0; i < no_blocks(); ++i)
    {
        M.block(i * m_block_rows,
                i * m_block_cols,
                m_block_rows,
                m_block_cols) = block(i);
    }

    return M;
}

BdMatrix BdMatrix::operator+ (const BdMatrix& rhs) const
{
    if (!dimensions_match(rhs))
        throw std::logic_error("Dimensions mismatch");

    BdMatrix result(no_blocks(), m_block_rows, m_block_cols);

    for (int i = 0; i < no_blocks(); ++i) {
        result.block(i) = block(i) + rhs.block(i);
    }

    return result;
}

BdMatrix BdMatrix::operator* (const BdMatrix& rhs) const
{
    if (no_blocks() != rhs.no_blocks() ||
            m_block_cols != rhs.m_block_rows)
        throw std::logic_error("Dimensions mismatch");

    BdMatrix result(no_blocks(), m_block_rows, rhs.m_block_cols);

    for (int i = 0; i < no_blocks(); ++i) {
        result.block(i) = block(i) * rhs.block(i);
    }

    return result;
}

BdMatrix operator* (complex<double> scalar, const BdMatrix& mat)
{
    BdMatrix result(mat.no_blocks(), mat.m_block_rows, mat.m_block_cols);

    for (int i = 0; i < mat.no_blocks(); ++i) {
        result.block(i) = scalar * mat.block(i);
    }

    return result;
}

bool BdMatrix::dimensions_match(const BdMatrix& rhs) const
{
    return (m_block_rows == rhs.m_block_rows)
            && (m_block_cols == rhs.m_block_cols)
            && (no_blocks() == rhs.no_blocks());
}

BdMatrix BdMatrix::inverse() const
{
    assert(m_block_rows == m_block_rows);
    BdMatrix result(no_blocks(), m_block_rows, m_block_cols);

    for (int i = 0; i < no_blocks(); ++i) {

        FullPivLU<MatrixXcd> lu(block(i));
        if (!lu.isInvertible()) {
            stringstream msg;
            msg << "Matrix is not invertible "
                << "(up to machine precision). "
                << "Failed to invert "
                << m_block_rows << "x" << m_block_cols
                << " block "
                << i << "." << endl
                << "Its numerical rank is "
                << lu.rank()
                << "." << endl;

            throw runtime_error(msg.str());
        }

        result.block(i) = lu.inverse();
    }

    return result;
}

BdMatrix BdMatrix::adjoint() const
{
    BdMatrix result(no_blocks(), m_block_cols, m_block_rows);

    for (int i = 0; i < no_blocks(); ++i) {
        result.block(i) = block(i).adjoint();
    }

    return result;
}

double BdMatrix::squaredNorm() const
{
    double sq_norm = 0;
    for (int i = 0; i < no_blocks(); ++i) {
        sq_norm += block(i).squaredNorm();
    }

    return sq_norm;
}

double BdMatrix::norm() const
{
    return sqrt(squaredNorm());
}

double BdMatrix::spectral_radius() const
{
    vector<double> radii(no_blocks());
    for (int i = 0; i < no_blocks(); ++i) {
        radii[i] = abs(eigenvalue_max_magnitude(block(i)));
    }

    return *std::max_element(radii.begin(), radii.end());
}

double BdMatrix::spectral_norm() const
{
    vector<double> norms(no_blocks());
    for (int i = 0; i < no_blocks(); ++i) {
        // compute the singular values
        norms[i] = sqrt(abs(eigenvalue_max_magnitude(
                    block(i).adjoint() * block(i))));
    }

    return *std::max_element(norms.begin(), norms.end());
}

VectorXcd BdMatrix::eigenvalues() const
{
    VectorXcd result(rows());
    int p = 0;
    for (int i = 0; i < no_blocks(); ++i) {
        // compute the eigenvalues of the block
        VectorXcd block_ews = lfa::eigenvalues(block(i));

        // add them to the vector of all eigenvalues
        for (int k = 0; k < block_rows(); ++k) {
            result(p) = block_ews(k);
            ++p;
        }
    }

    return result;
}

ostream& operator<< (ostream& os, const BdMatrix& m)
{
    for (int i = 0; i < m.no_blocks(); ++i) {
        os << "Block #" << i << endl;
        os << m.block(i) << endl;
    }
    return os;
}


}

