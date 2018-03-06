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

#ifndef LFA_SYMBOL_H
#define LFA_SYMBOL_H

#include "Common.h"
#include "BdMatrix.h"
#include "FrequencyIndices.h"
#include "HarmonicClusters.h"
#include "ClusterSymbol.h"

#include "CartIterator.h"
#include "Grid.h"
#include "SamplingProperties.h"
#include "NdArray.h"

namespace lfa {

/** Storage for the sampling of a block symbol. */
class Symbol {
    public:
        /** Create a symbol initialized with zeros. */
        Symbol(HarmonicClusters output_clusters = HarmonicClusters(),
               HarmonicClusters input_clusters = HarmonicClusters());

        static Symbol Identity(HarmonicClusters row_clusters,
                               HarmonicClusters col_clusters);
        static Symbol Identity(Grid grid, SamplingProperties conf);

        static Symbol Zero(HarmonicClusters row_clusters,
                           HarmonicClusters col_clusters);
        static Symbol Zero(Grid, SamplingProperties conf);

        BdMatrix& matrix() { return m_store; }

        Symbol addCompatible(const Symbol& other) const;
        Symbol operator+ (const Symbol& other) const;

        Symbol mulCompatible(const Symbol& other) const;
        Symbol operator* (const Symbol& other) const;

        friend Symbol operator* (complex<double> lhs, const Symbol& rhs);
        Symbol operator- (const Symbol& other) const {
            return (*this) + (-1.0 * other);
        }

        Symbol expand(ArrayFi factor) const;

        Symbol inverse() const;
        Symbol adjoint() const;

        complex<double>& ref(ArrayFi base,
                             ArrayFi cluster_row,
                             ArrayFi cluster_col);
        complex<double> ref(ArrayFi base,
                            ArrayFi cluster_row,
                            ArrayFi cluster_col) const {
            return const_cast<Symbol*>(this)->ref(base,
                                                  cluster_row,
                                                  cluster_col);
        }

        NdRange baseIndices() const { return m_output_clusters.baseIndices(); }

        friend class SymbolClusterRef;

        const HarmonicClusters& outputClusters() const { return m_output_clusters; }
        const HarmonicClusters& inputClusters() const { return m_input_clusters; }

        void setCluster(ArrayFi base, const ClusterSymbol& sym);
        ClusterSymbol getCluster(ArrayFi base) const;

        MatrixXcd full() const { return m_store.full(); }
        MatrixXcd fullCluster(ArrayFi b) const {
            return m_store.block(baseIndices().indexOf(b));
        }

        double norm() const;
        double squaredNorm() const;

        NdArray<double> row_norms() const;
        NdArray<double> col_norms() const;

        VectorXcd row_norms_1d() const;
        VectorXcd col_norms_1d() const;
        MatrixXcd row_norms_2d() const;
        MatrixXcd col_norms_2d() const;
        double spectral_radius() const;
        double spectral_norm() const;

        VectorXcd eigenvalues() const;

        int dimension() const { return m_output_clusters.dimension(); }

        class iterator;
        iterator begin();
        iterator end();
    private:
        HarmonicClusters m_output_clusters;

        HarmonicClusters m_input_clusters;
        BdMatrix m_store;
};

ostream& operator<<(ostream& os, Symbol sym);

class SymbolClusterRef {
    public:
        SymbolClusterRef(Symbol& symbol, ArrayFi base);

        complex<double>& operator() (ArrayFi cluster_row, ArrayFi cluster_col);

        void rowFrequency(ArrayFd step_size, ArrayFd base_freq, ArrayFi cluster_index);
    private:
        Symbol& m_symbol;
        ArrayFi m_base;
        int m_diag_index;
};

class Symbol::iterator : public std::iterator<std::input_iterator_tag,
                                              complex<double> >
{
    public:
        typedef CartIterator<NdRange::iterator, NdRange::iterator> InnerIterator;
        typedef CartIterator<InnerIterator, NdRange::iterator> OuterIterator;


        iterator(Symbol* sym, bool out_of_range = false)
        {
            m_sym = sym;

            NdRange c_grid = m_sym->inputClusters().clusterIndices();
            NdRange r_grid = m_sym->outputClusters().clusterIndices();

            InnerIterator inner_begin(c_grid.begin(), c_grid.end(),
                                      r_grid.begin(), r_grid.end());
            InnerIterator inner_end(  c_grid.begin(), c_grid.end(),
                                      r_grid.begin(), r_grid.end(), true);


            NdRange b_grid = m_sym->inputClusters().baseIndices();
            m_iter = OuterIterator(inner_begin, inner_end,
                                   b_grid.begin(), b_grid.end(), out_of_range);

        }

        bool operator!= (const iterator& other) {
            return m_iter != other.m_iter;
        }

        iterator& operator++ () {
            ++m_iter;

            return *this;
        }

        ArrayFi base() {
            return *(m_iter->second);
        }

        ArrayFi row() {
            return *(m_iter->first->second);
        }

        ArrayFi col() {
            return *(m_iter->first->first);
        }

        complex<double>& operator* () {
            return m_sym->ref( base(), row(), col() );
        }

    private:
        OuterIterator m_iter;

        Symbol* m_sym;
};

inline Symbol::iterator Symbol::begin() { return Symbol::iterator(this); }
inline Symbol::iterator Symbol::end() { return Symbol::iterator(this, true); }


}

#endif
