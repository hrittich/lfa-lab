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

#include "HarmonicClusters.h"
#include "MathUtil.h"

namespace lfa {

  HarmonicClusters::HarmonicClusters(ArrayFi base_shape,
      ArrayFi cluster_shape)
    : m_base_shape(base_shape),
    m_cluster_shape(cluster_shape)
  {

  }

  bool HarmonicClusters::operator== (const HarmonicClusters& other) const
  {
    return (m_cluster_shape == other.m_cluster_shape).all() &&
      (m_base_shape == other.m_base_shape).all();
  }

  ArrayFi HarmonicClusters::shape() const
  {
    return (m_base_shape * m_cluster_shape);
  }

  HarmonicClusters HarmonicClusters::mergeCluster(ArrayFi factor) const
  {
    if ((m_base_shape.binaryExpr(factor, std::modulus<int>())
          != ArrayFi::Zero(dimension())).all()) {
      throw logic_error("Base shape not sufficient.");
    }

    return HarmonicClusters(m_base_shape / factor, m_cluster_shape * factor);
  }


  int HarmonicClusters::clusterSize() const
  {
    return m_cluster_shape.prod();
  }

  int HarmonicClusters::baseSize() const
  {
    return m_base_shape.prod();
  }

  int HarmonicClusters::size() const
  {
    return baseSize() * clusterSize();
  }

  bool HarmonicClusters::isCompatibleTo(const HarmonicClusters& other) const
  {
    return (m_base_shape == other.m_base_shape).all();
  }

  ArrayFi HarmonicClusters::globalIndex(ArrayFi base_index, ArrayFi cluster_index) const
  {
    return (base_index + m_base_shape * cluster_index);
  }

  ArrayFi HarmonicClusters::baseIndex(ArrayFi global_index) const
  {
    return global_index.binaryExpr(m_base_shape, std::modulus<int>());
  }

  ArrayFi HarmonicClusters::clusterIndex(ArrayFi global_index) const
  {
    return global_index / m_base_shape;
  }

  void HarmonicClusters::convert(ArrayFi& result_base_index,
      ArrayFi& result_cluster_index,
      HarmonicClusters result_clusters,
      ArrayFi base_index,
      ArrayFi cluster_index) const
  {
    ArrayFi global = globalIndex(base_index, cluster_index);
    result_base_index = result_clusters.baseIndex(global);
    result_cluster_index = result_clusters.clusterIndex(global);
  }


  HarmonicClusters HarmonicClusters::minContainer(const HarmonicClusters& other) const
  {
    ArrayFi l = lcm(this->clusterIndices().shape(),
        other.clusterIndices().shape());

    if ( (shape().binaryExpr(l, std::modulus<int>())
          != ArrayFi::Zero(dimension())).any() )
      throw logic_error("Container size exhausted");

    return HarmonicClusters(shape() / l, l);
  }

  ArrayFi HarmonicClusters::expansionFactor(HarmonicClusters& expanded) const
  {
    if ((expanded.m_cluster_shape.binaryExpr(m_cluster_shape, std::modulus<int>())
          != ArrayFi::Zero(dimension())).any()) {
      throw logic_error("Cannot expand");
    }

    return expanded.m_cluster_shape / m_cluster_shape;
  }


}

