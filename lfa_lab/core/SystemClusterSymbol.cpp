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

#include "SystemClusterSymbol.h"
#include <Eigen/Dense>

namespace lfa {

SystemClusterSymbol::SystemClusterSymbol(int rows,
                                         int cols,
                                         ArrayFi row_cluster_shape,
                                         ArrayFi col_cluster_shape)
  : m_rows(rows),
    m_cols(cols),
    m_row_indices(row_cluster_shape),
    m_col_indices(col_cluster_shape)
{
    m_store.resize(rows * m_row_indices.size(),
                   cols * m_col_indices.size());
}

void SystemClusterSymbol::setCluster(int si, int sj, const ClusterSymbol& cluster)
{
    assert( (m_row_indices.shape() == cluster.rowShape()).all() );
    assert( (m_col_indices.shape() == cluster.colShape()).all() );

    for (NdRange::iterator pi = m_row_indices.begin();
            pi != m_row_indices.end(); ++pi)
    {
        for (NdRange::iterator pj = m_col_indices.begin();
                pj != m_col_indices.end(); ++pj)
        {
            (*this)(si, sj, *pi, *pj) = cluster(*pi, *pj);
        }
    }
}

ClusterSymbol SystemClusterSymbol::getCluster(int si, int sj)
{
    ClusterSymbol cluster(m_row_indices.shape(),
                          m_col_indices.shape());

    for (NdRange::iterator pi = m_row_indices.begin();
            pi != m_row_indices.end(); ++pi)
    {
        for (NdRange::iterator pj = m_col_indices.begin();
                pj != m_col_indices.end(); ++pj)
        {
            cluster(*pi, *pj) = (*this)(si, sj, *pi, *pj);
        }
    }

    return cluster;
}

SystemClusterSymbol SystemClusterSymbol::inverse() const
{
    SystemClusterSymbol result(m_cols,
                               m_rows,
                               m_col_indices.shape(),
                               m_row_indices.shape());

    FullPivLU<MatrixXcd> lu(m_store);
    if (!lu.isInvertible()) {
        throw runtime_error("Matrix is not invertible "
                "(up to machine precision)");
    }
    result.m_store = lu.inverse();

    return result;
}

complex<double>& SystemClusterSymbol::operator() (int si,
                                                  int sj,
                                                  const ArrayFi& ci,
                                                  const ArrayFi& cj)
{
    assert(0 <= si && si < m_rows);
    assert(0 <= sj && sj < m_cols);

    const int m_row_stride = m_row_indices.size();
    const int m_col_stride = m_col_indices.size();

    int i = m_row_stride * si + m_row_indices.indexOf(ci);
    int j = m_col_stride * sj + m_col_indices.indexOf(cj);

    return m_store(i, j);
}



}

