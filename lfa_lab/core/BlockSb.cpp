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

#include "BlockSb.h"
#include "DiscreteDomain.h"
#include "MathUtil.h"

namespace lfa {

BlockSb::BlockSb(Grid grid, ArrayFi period)
  : m_grid(grid), m_period(period)
{
}

void BlockSb::scalarSymbols(NdArray<Symbol> scalars)
{
    if ((scalars.shape() != m_period).any()) {
        throw logic_error("Invalid number of scalars.");
    }

    m_scalars = scalars;
}

FoProperties BlockSb::properties()
{
    SplitFrequencyDomain domain(m_grid, m_period);
    FoProperties prop(domain, domain);

    return prop;
}

SplitFrequencyDomain BlockSb::scalar_domain()
{
    return SplitFrequencyDomain(m_grid, ArrayFi::Ones(dimension()));
}

FoProperties BlockSb::scalar_properties()
{
    FoProperties prop(scalar_domain(), scalar_domain());

    return prop;
}

Symbol BlockSb::generate(const SamplingProperties& conf)
{
    DiscreteDomain scalar_dd(scalar_domain(), conf);

    // we need to ensure that all symbols habe the same properties
    NdRange indices = m_scalars.indices();
    for (NdRange::iterator p = indices.begin(); p != indices.end(); ++p) {
        if (m_scalars(*p).inputClusters() != scalar_dd.harmonics()
                || m_scalars(*p).outputClusters() != scalar_dd.harmonics())
        {
            throw logic_error("Symbol has the wrong coupling.");
        }
    }


    // do the computation
    DiscreteDomain domain(properties().output(), conf);

    ArrayFi zero = ArrayFi::Zero(dimension());
    HarmonicClusters clusters = domain.harmonics();

    NdRange base_grid = clusters.baseIndices();
    NdRange cluster_grid = clusters.clusterIndices();

    Symbol F(clusters, clusters);
    Symbol G(clusters, clusters);

    for (NdRange::iterator b = base_grid.begin(); b != base_grid.end(); ++b)
    {
        for (NdRange::iterator i = cluster_grid.begin();
                i != cluster_grid.end(); ++i)
        {
            for (NdRange::iterator j = cluster_grid.begin();
                    j != cluster_grid.end(); ++j)
            {
                ArrayFi gj = clusters.globalIndex(*b, *j);

                double arg = 2 * pi *
                             (i->cast<double>() *
                              j->cast<double>() /
                              clusters.clusterShape().cast<double>())
                             .sum();
                F.ref(*b, *i, *j) = exp(complex<double>(0, arg));
                G.ref(*b, *i, *j) = m_scalars(*i).ref(gj,zero,zero)
                                    * F.ref(*b, *i, *j);

            }
        }
    }

    return (1.0 / clusters.clusterSize()) * F.adjoint() * G;
}



}

