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

#include "ConstantSb.h"
#include "DiscreteDomain.h"

namespace lfa {

ConstantSb::ConstantSb()
{

}

ConstantSb::ConstantSb(ClusterSymbol symbol, Grid output_grid, Grid input_grid)
  : m_symbol(symbol),
    m_output_grid(output_grid),
    m_input_grid(input_grid)
{

}

FoProperties ConstantSb::properties()
{
    return FoProperties(
        SplitFrequencyDomain(m_output_grid, m_symbol.outputShape()),
        SplitFrequencyDomain(m_input_grid, m_symbol.inputShape()));
}

Symbol ConstantSb::generate(const SamplingProperties& conf)
{
    Symbol result(
            DiscreteDomain(properties().output(), conf).harmonics(),
            DiscreteDomain(properties().input(), conf).harmonics());

    NdRange indices = result.baseIndices();
    for (NdRange::iterator b = indices.begin(); b != indices.end(); ++b)
    {
        result.setCluster(*b, m_symbol);
    }

    return result;
}

ClusterSymbol flat_restriction_cluster_symbol(Grid output_grid, Grid input_grid)
{
    ArrayFi factor = input_grid.coarsening_factor(output_grid);
    int d = factor.rows();

    ClusterSymbol sym(ArrayFi::Ones(d), factor);

    NdRange indices(factor);
    double s = sqrt(1.0 / indices.size());
    for (NdRange::iterator p = indices.begin();
            p != indices.end(); ++p)
    {
        sym(ArrayFi::Zero(d), *p) = s;
    }

    return sym;
}

ConstantSb flat_restriction_sb(Grid output_grid, Grid input_grid)
{
    return ConstantSb(
            flat_restriction_cluster_symbol(output_grid, input_grid),
            output_grid, input_grid);
}

ClusterSymbol flat_interpolation_cluster_symbol(Grid output_grid, Grid input_grid)
{
    ArrayFi factor = output_grid.coarsening_factor(input_grid);
    int d = factor.rows();

    ClusterSymbol sym(factor, ArrayFi::Ones(d));

    NdRange indices(factor);
    double s = sqrt(1.0 / indices.size());

    for (NdRange::iterator p = indices.begin();
            p != indices.end(); ++p)
    {
        sym(*p, ArrayFi::Zero(d)) = s;
    }

    return sym;
}

ConstantSb flat_interpolation_sb(Grid output_grid, Grid input_grid)
{
    return ConstantSb(
            flat_interpolation_cluster_symbol(output_grid, input_grid),
            output_grid, input_grid);
}

ClusterSymbol zero_restriction_cluster_symbol(Grid output_grid, Grid input_grid)
{
    ArrayFi factor = input_grid.coarsening_factor(output_grid);
    int d = factor.rows();

    ClusterSymbol sym(ArrayFi::Ones(d), factor);

    NdRange indices(factor);
    for (NdRange::iterator p = indices.begin();
            p != indices.end(); ++p)
    {
        sym(ArrayFi::Zero(d), *p) = 0;
    }

    return sym;
}

ConstantSb zero_restriction_sb(Grid output_grid, Grid input_grid)
{
    return ConstantSb(
            zero_restriction_cluster_symbol(output_grid, input_grid),
            output_grid, input_grid);
}

ClusterSymbol zero_interpolation_cluster_symbol(Grid output_grid, Grid input_grid)
{
    ArrayFi factor = output_grid.coarsening_factor(input_grid);
    int d = factor.rows();

    ClusterSymbol sym(factor, ArrayFi::Ones(d));

    NdRange indices(factor);
    for (NdRange::iterator p = indices.begin();
            p != indices.end(); ++p)
    {
        sym(*p, ArrayFi::Zero(d)) = 0;
    }

    return sym;
}

ConstantSb zero_interpolation_sb(Grid output_grid, Grid input_grid)
{
    return ConstantSb(
            zero_interpolation_cluster_symbol(output_grid, input_grid),
            output_grid, input_grid);
}



}


