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

#ifndef LFA_SYSTEM_CLUSTER_SYMBOL_H
#define LFA_SYSTEM_CLUSTER_SYMBOL_H

#include "Common.h"
#include "ClusterSymbol.h"

namespace lfa {

  class SystemClusterSymbol {
    public:
      SystemClusterSymbol(int rows,
                          int cols,
                          ArrayFi row_cluster_shape,
                          ArrayFi col_cluster_shape);

      void setCluster(int si, int sj, const ClusterSymbol& cluster);
      ClusterSymbol getCluster(int si, int sj);

      SystemClusterSymbol inverse() const;

      MatrixXcd& matrix() { return m_store; }
    private:
      complex<double>& operator() (int si,
                                   int sj,
                                   const ArrayFi& ci,
                                   const ArrayFi& cj);

      int m_rows;
      int m_cols;
      NdRange m_row_indices;
      NdRange m_col_indices;
      MatrixXcd m_store;
  };

}

#endif