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

#ifndef LFA_CLUSTER_SYMBOL_H
#define LFA_CLUSTER_SYMBOL_H

#include "Common.h"
#include "NdRange.h"

namespace lfa {

  /** The matrix representing the coupling of a set of harmonic frequencies.
   */
  class ClusterSymbol {
    public:
      ClusterSymbol(ArrayFi output_shape = ArrayFi::Zero(0),
                    ArrayFi input_shape = ArrayFi::Zero(0));

      complex<double>& operator() (ArrayFi row, ArrayFi col);
      complex<double> operator() (ArrayFi row, ArrayFi col) const {
        return (*const_cast<ClusterSymbol*>(this))(row, col);
      }

      MatrixXcd toMatrix() const { return m_store; }
      void setMatrix(MatrixXcd m);

      ArrayFi outputShape() const { return m_output_indices.shape(); }
      ArrayFi inputShape() const { return m_input_indices.shape(); }

      ArrayFi rowShape() const { return rowIndices().shape(); }
      ArrayFi colShape() const { return colIndices().shape(); }

      NdRange rowIndices() const { return m_output_indices; }
      NdRange colIndices() const { return m_input_indices; }
    private:
      NdRange m_output_indices;
      NdRange m_input_indices;
      MatrixXcd m_store;
  };

};

#endif
