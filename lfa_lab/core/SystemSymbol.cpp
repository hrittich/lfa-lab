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

#include "SystemSymbol.h"
#include "SystemClusterSymbol.h"
#include "ExEigenSolver.h"

namespace lfa {

SystemSymbol::SystemSymbol(int rows,
                           int cols,
                           HarmonicClusters output_clusters,
                           HarmonicClusters input_clusters)
  : m_output_clusters(output_clusters),
    m_input_clusters(input_clusters)
{
  resize(rows, cols);
}

SystemSymbol SystemSymbol::Identity(int rows,
                                    int cols,
                                    HarmonicClusters output_clusters,
                                    HarmonicClusters input_clusters)
{
  SystemSymbol result(rows, cols, output_clusters, input_clusters);

  Symbol I = Symbol::Identity(output_clusters, input_clusters);
  Symbol Z = Symbol::Zero(output_clusters, input_clusters);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (i == j) {
        result(i, j) = I;
      } else {
        result(i, j) = Z;
      }
    }
  }

  return result;
}

void SystemSymbol::resize(int rows, int cols)
{
  m_rows = rows;
  m_cols = cols;

  vector<Symbol> default_row(cols);
  std::fill(default_row.begin(),
            default_row.end(),
            Symbol(m_output_clusters, m_input_clusters));
  m_store.resize(rows);
  std::fill(m_store.begin(), m_store.end(), default_row);
}

SystemSymbol SystemSymbol::operator* (const SystemSymbol& other) const
{
  assert( cols() == other.rows() );

  #ifndef NDEBUG
  ensureConsistency();
  other.ensureConsistency();
  #endif

  // find the common cluster for the output of the first and the input of the
  // second
  HarmonicClusters common
    = m_input_clusters.minContainer(other.m_output_clusters);

  ArrayFi self_factor = m_input_clusters.expansionFactor(common);
  ArrayFi other_factor = other.m_output_clusters.expansionFactor(common);

  HarmonicClusters output_clusters
    = m_output_clusters.mergeCluster(self_factor);
  HarmonicClusters input_clusters
    = other.m_input_clusters.mergeCluster(other_factor);

  SystemSymbol result(m_rows,
                      other.m_cols,
                      output_clusters,
                      input_clusters);

  for (int i = 0; i < result.rows(); ++i) {
    for (int j = 0; j < result.cols(); ++j) {
      Symbol aux = Symbol::Zero(output_clusters, input_clusters);

      for (int k = 0; k < cols(); ++k) {
        aux = aux + (*this)(i, k) * other(k, j);
      }

      result(i,j) = aux;
    }
  }

  return result;
}

SystemSymbol SystemSymbol::operator+ (const SystemSymbol& other) const
{
  HarmonicClusters common_input =
    m_input_clusters.minContainer(other.m_input_clusters);

  ArrayFi factor = m_input_clusters.expansionFactor(common_input);
  HarmonicClusters common_output = m_output_clusters.mergeCluster(factor);


  SystemSymbol result(m_rows, m_cols, common_output, common_input);

  for (int i = 0; i < m_rows; ++i) {
    for (int j = 0; j < m_cols; ++j) {
      result(i,j) = (*this)(i,j) + other(i,j);
    }
  }

  return result;
}

SystemSymbol operator* (double scalar, const SystemSymbol& other)
{
  SystemSymbol result(other.m_rows,
                      other.m_cols,
                      other.m_output_clusters,
                      other.m_input_clusters);

  for (int i = 0; i < other.m_rows; ++i) {
    for (int j = 0; j < other.m_cols; ++j) {
      result(i,j) = scalar * other(i,j);
    }
  }
  return result;
}

SystemClusterSymbol SystemSymbol::at(ArrayFi base_index) const
{
  ensureConsistency();

  SystemClusterSymbol aux(m_rows,
        m_cols,
        m_output_clusters.clusterShape(),
        m_input_clusters.clusterShape());

  for (int si = 0; si < m_rows; ++si) {
    for (int sj = 0; sj < m_cols; ++sj) {
      aux.setCluster(si, sj, (*this)(si, sj).getCluster(base_index));
    }
  }

  return aux;
}

SystemSymbol SystemSymbol::inverse() const
{
  SystemSymbol result(m_cols, m_rows, m_input_clusters, m_output_clusters);

  NdRange bases = baseIndices();
  for (NdRange::iterator b = bases.begin(); b != bases.end(); ++b)
  {
    SystemClusterSymbol aux = at(*b);

    SystemClusterSymbol aux_inv = aux.inverse();

    // Assert that we really computed the inverse
    MatrixXcd expected_id = aux_inv.matrix() * aux.matrix();
    if ( (expected_id -
          MatrixXcd::Identity(aux.matrix().rows(),
            aux.matrix().cols())).norm()
        >= 1e-12 )
    {
      throw runtime_error("Inversion failed");
    }

    for (int si = 0; si < m_rows; ++si) {
      for (int sj = 0; sj < m_cols; ++sj) {
        result(si, sj).setCluster(*b, aux_inv.getCluster(si, sj));
      }
    }
  }

  return result;
}

double SystemSymbol::spectral_radius() const
{
  vector<double> max_evs;

  NdRange bases = baseIndices();
  for (NdRange::iterator b = bases.begin(); b != bases.end(); ++b)
  {
    SystemClusterSymbol aux = at(*b);

    max_evs.push_back(
      abs(eigenvalue_max_magnitude(aux.matrix())));
  }

  return *max_element(max_evs.begin(), max_evs.end());
}

double SystemSymbol::spectral_norm() const
{
  vector<double> max_evs;

  NdRange bases = baseIndices();
  for (NdRange::iterator b = bases.begin(); b != bases.end(); ++b)
  {
    SystemClusterSymbol aux = at(*b);

    max_evs.push_back(
      abs(eigenvalue_max_magnitude(aux.matrix().adjoint()
                                   * aux.matrix())));
  }

  return *max_element(max_evs.begin(), max_evs.end());
}

double SystemSymbol::system_norm() const
{
  double norm_sq = 0;

  for (int i = 0; i < rows(); ++i) {
    for (int j = 0; j < cols(); ++j) {
      norm_sq += (*this)(i,j).squaredNorm();
    }
  }
  return sqrt(norm_sq);
}

void SystemSymbol::ensureConsistency() const
{
  HarmonicClusters output_clusters = (*this)(0,0).outputClusters();
  HarmonicClusters input_clusters = (*this)(0,0).inputClusters();

  // ensure that the input and output clusters are all identical
  for (int i=0; i < rows(); ++i) {
    for (int j=0; j < cols(); ++j) {
      if ((*this)(i,j).inputClusters() != input_clusters ||
          (*this)(i,j).outputClusters() != output_clusters)
      {
        throw std::runtime_error("Symbol need to have the same clusters");
      }
    }
  }
}

SystemSymbol combine_symbols_into_system(
  const SystemSymbolProperties& properties,
  const MatrixContainer<Symbol>& symbols)
{
  if (symbols.rows() <= 0 || symbols.cols() <= 0) {
    throw std::logic_error("Rows or columns of symbol matrix are zero, which "
                           "is not allowed");
  }

  HarmonicClusters output_clusters
    = make_harmonic_cluster(symbols(0,0).outputClusters().shape(),
                            properties.element_properties().output());

  HarmonicClusters input_clusters
    = make_harmonic_cluster(symbols(0,0).inputClusters().shape(),
                            properties.element_properties().input());

  SystemSymbol symbol(symbols.rows(),
                      symbols.cols(),
                      output_clusters,
                      input_clusters);

  for (int i=0; i < symbols.rows(); ++i) {
    for (int j=0; j < symbols.cols(); ++j) {
      ArrayFi factor = input_clusters.clusterShape()
                        / symbols(i,j).inputClusters().clusterShape();
      symbol(i,j) = symbols(i,j).expand(factor);
    }
  }

  symbol.ensureConsistency();

  return symbol;
}



}

