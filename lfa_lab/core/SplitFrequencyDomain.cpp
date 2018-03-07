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

#include "SplitFrequencyDomain.h"
#include "MathUtil.h"

namespace lfa {

  SplitFrequencyDomain::SplitFrequencyDomain(Grid grid, ArrayFi cluster_shape)
    : m_grid(grid), m_cluster_shape(cluster_shape)
  {
    if (m_cluster_shape.rows() == 0) {
      m_cluster_shape = ArrayFi::Ones(m_grid.dimension());
    }
  }

  HarmonicClusters SplitFrequencyDomain::harmonics(ArrayFi resolution) const
  {
    ArrayFi zero = ArrayFi::Zero(resolution.rows());

    if ((resolution.binaryExpr(m_cluster_shape, std::modulus<int>())
          != zero).any())
    {
      stringstream msg;
      msg << "You requested a resolution of " << ListFmt(resolution)
        << " however it has to be a multiple of " << ListFmt(m_cluster_shape)
        << ".";
      throw logic_error(msg.str());
    }

    return HarmonicClusters(
        resolution / m_cluster_shape,
        m_cluster_shape);
  }

  ArrayFi SplitFrequencyDomain::resolution(ArrayFi finest_resolution) const
  {
    assert( (finest_resolution.binaryExpr(m_grid.spacing(), std::modulus<int>())
          == ArrayFi::Zero(dimension())).all() );

    return  finest_resolution / m_grid.spacing();

  }

  bool SplitFrequencyDomain::isValid()
  {
    return
      ((m_cluster_shape > ArrayFi::Zero(m_cluster_shape.rows())).all()
       && m_cluster_shape.rows() > 0);
  }

  bool SplitFrequencyDomain::isCompatibleTo(const SplitFrequencyDomain& other)
  {
    return (m_grid.spacing() * m_cluster_shape ==
        other.m_grid.spacing() * other.m_cluster_shape).all();
  }

  SplitFrequencyDomain SplitFrequencyDomain::expand(ArrayFi factor) const
  {
    return SplitFrequencyDomain(m_grid, m_cluster_shape * factor);
  }

  ArrayFi SplitFrequencyDomain::lcc(const SplitFrequencyDomain& other) const
  {
    return lcm(m_cluster_shape, other.m_cluster_shape);
  }

  ostream& operator<< (ostream& os, const SplitFrequencyDomain& domain) {
    os << "SplitFrequencyDomain("
      << "grid = " << domain.grid() << ", "
      << "cluster_shape = " << ListFmt(domain.clusterShape()) << ")";

    return os;
  }



}

