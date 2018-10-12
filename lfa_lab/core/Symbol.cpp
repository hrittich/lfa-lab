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

#include "Symbol.h"

#include "FrequencyIndices.h"
#include "MathUtil.h"
#include "SplitFrequencyDomain.h"
#include "DiscreteDomain.h"

namespace lfa {

    Symbol::Symbol(HarmonicClusters output_clusters, HarmonicClusters input_clusters)
        : m_output_clusters(output_clusters),
          m_input_clusters(input_clusters)
    {
        // the input and output harmonics have to have the same base indices
        if (!m_input_clusters.isCompatibleTo(m_output_clusters))
            throw logic_error("Input and output modes are incompatibel.");

        m_store.resize(
            output_clusters.baseIndices().size(),
            output_clusters.clusterIndices().size(),
            input_clusters.clusterIndices().size());
        m_store.setZero();
    }

    Symbol Symbol::Identity(HarmonicClusters output_clusters,
                            HarmonicClusters input_clusters)
    {
        Symbol result(output_clusters, input_clusters);
        for (Symbol::iterator p = result.begin(); p != result.end(); ++p)
        {
            if ( (p.row() == p.col()).all() ) {
                *p = 1;
            } else {
                *p = 0;
            }
        }
        return result;
    }

    Symbol Symbol::Identity(Grid grid, SamplingProperties conf)
    {
        SplitFrequencyDomain cont_domain(grid);
        DiscreteDomain domain(cont_domain, conf);
        return Identity(domain.harmonics(), domain.harmonics());
    }

    Symbol Symbol::Zero(HarmonicClusters output_clusters,
                        HarmonicClusters input_clusters)
    {
        return Symbol(output_clusters, input_clusters);
    }

    Symbol Symbol::Zero(Grid grid, SamplingProperties conf)
    {
        SplitFrequencyDomain cont_domain(grid);
        DiscreteDomain domain(cont_domain, conf);
        return Zero(domain.harmonics(), domain.harmonics());
    }


    Symbol Symbol::addCompatible(const Symbol& other) const
    {
        if ( m_input_clusters != other.m_input_clusters
             || m_output_clusters != other.m_output_clusters)
            throw logic_error("Symbols do not correspond to the same harmonics");

        Symbol result(m_output_clusters, m_input_clusters);
        result.m_store = m_store + other.m_store;

        return result;
    }

    Symbol Symbol::operator+ (const Symbol& other) const
    {
        HarmonicClusters common_input =
            m_input_clusters.minContainer(other.m_input_clusters);

        Symbol first = this->expand(m_input_clusters.expansionFactor(common_input));
        Symbol second = other.expand(other.m_input_clusters.expansionFactor(common_input));

        return first.addCompatible(second);
    }

    Symbol Symbol::mulCompatible(const Symbol& other) const
    {
        if ( m_input_clusters != other.m_output_clusters )
            throw logic_error("Symbols not compatible");

        Symbol result(m_output_clusters, other.m_input_clusters);
        result.m_store = m_store * other.m_store;

        return result;
    }

    Symbol Symbol::operator* (const Symbol& other) const
    {
        HarmonicClusters common =
            m_input_clusters.minContainer(other.m_output_clusters);

        Symbol first = this->expand(m_input_clusters.expansionFactor(common));
        Symbol second = other.expand(other.m_output_clusters.expansionFactor(common));

        return first.mulCompatible(second);
    }


    Symbol operator* (complex<double> lhs, const Symbol& rhs)
    {
        Symbol result(rhs.m_output_clusters, rhs.m_input_clusters);
        result.m_store = lhs * rhs.m_store;

        return result;
    }

    Symbol Symbol::expand(ArrayFi factor) const
    {
        // the input and output harmonics have to have the same base indices
        if (!m_input_clusters.isCompatibleTo(m_output_clusters))
            throw logic_error("Input and output modes are incompatibel.");

        // check if we need to expand
        if (factor.isConstant(1)) {
            return *this;
        }

        Symbol result(m_output_clusters.mergeCluster(factor),
                      m_input_clusters.mergeCluster(factor));
        result.m_store.setZero();

        // run through the small blocks an distribute the entries
        NdRange base_idx = m_input_clusters.baseIndices();
        NdRange row_cluster_idx = m_output_clusters.clusterIndices();
        NdRange col_cluster_idx = m_input_clusters.clusterIndices();

        // for each diagonal block
        for (NdRange::iterator b = base_idx.begin();
             b != base_idx.end(); ++b)
        {
            // for all rows
            for (NdRange::iterator rc = row_cluster_idx.begin();
                 rc != row_cluster_idx.end(); ++rc)
            {
                // for all columns
                for (NdRange::iterator cc = col_cluster_idx.begin();
                     cc != col_cluster_idx.end(); ++cc)
                {
                    // convert to target coordinates
                    ArrayFi target_rb, target_rc;
                    m_output_clusters.convert(target_rb, target_rc,
                                              result.m_output_clusters, *b, *rc);

                    ArrayFi target_cb, target_cc;
                    m_input_clusters.convert(target_cb, target_cc,
                                             result.m_input_clusters, *b, *cc);

                    assert((target_rb == target_cb).all());

                    // copy value
                    result.ref(target_rb, target_rc, target_cc) = ref(*b, *rc, *cc);
                }
            }
        }

        return result;
    }

    Symbol Symbol::inverse() const
    {
        assert(m_output_clusters == m_input_clusters);
        Symbol result(m_output_clusters, m_input_clusters);
        result.m_store = m_store.inverse();
        return result;
    }

    Symbol Symbol::adjoint() const
    {
        Symbol result(m_input_clusters, m_output_clusters);
        result.m_store = m_store.adjoint();
        return result;
    }

    complex<double>& Symbol::ref(ArrayFi base, ArrayFi cluster_row, ArrayFi cluster_col)
    {
        // the input and output harmonics have to have the same base indices
        if (!m_input_clusters.isCompatibleTo(m_output_clusters))
            throw logic_error("Input and output modes are incompatibel.");

        SymbolClusterRef cluster_ref(*this, base);

        return cluster_ref(cluster_row, cluster_col);
    }

    void Symbol::setCluster(ArrayFi base, const ClusterSymbol& sym)
    {
        assert( (m_output_clusters.clusterShape() == sym.rowShape()).all() );
        assert( (m_input_clusters.clusterShape() == sym.colShape()).all() );

        int b = baseIndices().indexOf(base);
        m_store.set_block(b, sym.toMatrix());
    }

    ClusterSymbol Symbol::getCluster(ArrayFi base) const
    {
        int b = baseIndices().indexOf(base);

        ClusterSymbol result(m_output_clusters.clusterShape(),
                             m_input_clusters.clusterShape());
        result.setMatrix(m_store.block(b));

        return result;
    }

    double Symbol::norm() const {
        return m_store.norm();
    }

    double Symbol::squaredNorm() const {
        return m_store.squaredNorm();
    }

    NdArray<double> Symbol::row_norms() const
    {
        NdRange output_grid = outputClusters().globalIndices();

        NdArray<double> result_sq(output_grid.shape());
        for (NdRange::iterator i = output_grid.begin();
             i != output_grid.end(); ++i)
        {
            result_sq(*i) = 0;
        }

        for (iterator p = const_cast<Symbol*>(this)->begin();
             p != const_cast<Symbol*>(this)->end(); ++p)
        {
            ArrayFi i = outputClusters().globalIndex(p.base(), p.row());
            //ArrayFi j = inputClusters().globalIndex(p.base(), p.col());

            result_sq(i) += abs_sq(*p);
        }

        NdArray<double> result(output_grid.shape());
        for (NdRange::iterator i = output_grid.begin();
             i != output_grid.end(); ++i)
        {
            result(*i) = sqrt(result_sq(*i));
        }

        return result;
    }

    NdArray<double> Symbol::col_norms() const {
        NdRange input_grid = inputClusters().globalIndices();

        NdArray<double> result_sq(input_grid.shape());
        for (NdRange::iterator j = input_grid.begin();
             j != input_grid.end(); ++j)
        {
            result_sq(*j) = 0;
        }

        for (iterator p = const_cast<Symbol*>(this)->begin();
             p != const_cast<Symbol*>(this)->end(); ++p)
        {
            ArrayFi j = inputClusters().globalIndex(p.base(), p.col());
            result_sq(j) += abs_sq(*p);
        }

        NdArray<double> result(input_grid.shape());
        for (NdRange::iterator j = input_grid.begin();
             j != input_grid.end(); ++j)
        {
            result(*j) = sqrt(result_sq(*j));
        }

        return result;
    }

    MatrixXcd Symbol::row_norms_2d() const
    {
        return to_matrix(row_norms());
    }

    MatrixXcd Symbol::col_norms_2d() const
    {
        return to_matrix(col_norms());
    }

    VectorXcd Symbol::row_norms_1d() const
    {
        return to_vector(row_norms());
    }

    VectorXcd Symbol::col_norms_1d() const
    {
        return to_vector(col_norms());
    }

    double Symbol::spectral_radius() const
    {
        return m_store.spectral_radius();
    }

    double Symbol::spectral_norm() const
    {
        return m_store.spectral_norm();
    }

    VectorXcd Symbol::eigenvalues() const
    {
        return m_store.eigenvalues();
    }

    ostream& operator<<(ostream& os, Symbol sym) {
        os << sym.matrix();
        return os;
    }

    SymbolClusterRef::SymbolClusterRef(Symbol& symbol, ArrayFi base)
        : m_symbol(symbol),
          m_base(base)
    {
        // the input and output harmonics have to have the same base indices
        if (!symbol.m_input_clusters.isCompatibleTo(symbol.m_output_clusters))
            throw logic_error("Input and output modes are incompatibel.");

        m_diag_index = symbol.m_output_clusters.baseIndices().indexOf(base);
    }

    complex<double>& SymbolClusterRef::operator() (ArrayFi cluster_row, ArrayFi cluster_col)
    {
        int irc = m_symbol.m_output_clusters.clusterIndices().indexOf(cluster_row);
        int icc = m_symbol.m_input_clusters.clusterIndices().indexOf(cluster_col);

        return m_symbol.m_store(m_diag_index, irc, icc);
    }

}

