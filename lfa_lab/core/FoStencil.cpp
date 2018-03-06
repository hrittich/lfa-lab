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

#include "FoStencil.h"

#include "SplitFrequencyDomain.h"
#include "DiscreteDomain.h"

namespace lfa {

  FoStencil::FoStencil(const SparseStencil& stencil, Grid grid)
    : m_stencil(stencil),
    m_grid(grid)
  {

  }

  FoProperties FoStencil::properties()
  {
    SplitFrequencyDomain domain(
        m_grid,
        ArrayFi::Ones(m_grid.dimension()));

    return FoProperties(domain, domain);
  }

  Symbol FoStencil::generate(const SamplingProperties& conf)
  {
    SplitFrequencyDomain cont_domain(
        m_grid,
        ArrayFi::Ones(m_grid.dimension()));

    DiscreteDomain domain(cont_domain, conf);

    HarmonicClusters clusters = domain.harmonics();
    Symbol sym(clusters, clusters);

    NdRange bases = clusters.baseIndices();
    for (NdRange::iterator b = bases.begin(); b != bases.end(); ++b)
    {
      SymbolClusterRef cluster(sym, *b);
      fill(cluster, *b, domain);
    }

    return sym;
  }

  void FoStencil::fill(SymbolClusterRef cluster,
                       ArrayFi base_index,
                       const DiscreteDomain& domain)
  {
    ArrayFi zero = ArrayFi::Zero(dimension());
    VectorFd frequency = domain.frequency(base_index, zero);
    cluster(zero, zero) = symbolAt(frequency);
  }

  complex<double> FoStencil::symbolAt(VectorFd frequency)
  {
    complex<double> m = 0;
    for (SparseStencil::iterator it = m_stencil.begin();
        it != m_stencil.end(); ++it) {
      VectorFd pos = it->offset.cast<double>() * m_grid.step_size();

      m += it->value * exp(complex<double>(0, frequency.dot(pos) ));
    }

    return m;
  }

}

