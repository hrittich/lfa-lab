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

#include "FoProperties.h"

namespace lfa {

  FoProperties::FoProperties(
      SplitFrequencyDomain output_domain,
      SplitFrequencyDomain input_domain)
    : m_output_domain(output_domain),
    m_input_domain(input_domain)
  {

  }

  FoProperties FoProperties::operator+ (const FoProperties& other)
  {
    ArrayFi lcc = m_input_domain.lcc(other.m_input_domain);
    FoProperties this_ex = this->expandInputTo(lcc);
    FoProperties other_ex = other.expandInputTo(lcc);

    if (this_ex.m_output_domain != other_ex.m_output_domain
        || this_ex.m_input_domain != other_ex.m_input_domain) {
      throw logic_error("Incompatible addition");
    }

    return this_ex;
  }

  FoProperties FoProperties::operator* (const FoProperties& other)
  {
    // ToDo simplify this
    ArrayFi lcc = m_input_domain.lcc(other.m_output_domain);
    FoProperties this_ex = this->expandInputTo(lcc);
    FoProperties other_ex = other.expandOutputTo(lcc);

    if (this_ex.m_input_domain != other_ex.m_output_domain) {
      stringstream msg;
      msg << "Incompatible contraction dimension. "
        << "The following domains mismatch. " << endl
        << "First input:   " << this_ex.m_input_domain << endl
        << "Second output: " << other_ex.m_output_domain;

      throw logic_error(msg.str());
    }

    return FoProperties(this_ex.m_output_domain, other_ex.m_input_domain);
  }

  FoProperties FoProperties::inverse()
  {
    return FoProperties(m_input_domain, m_output_domain);
  }

  FoProperties FoProperties::adjoint()
  {
    return FoProperties(m_input_domain, m_output_domain);
  }

  ArrayFi FoProperties::adjustResolution(ArrayFi desired_resolution)
  {
    ArrayFi factor =
      m_output_domain.grid().spacing() * m_output_domain.clusterShape();

    ArrayFi rem =
      desired_resolution.binaryExpr(factor, std::modulus<int>());

    if ( (rem == ArrayFi::Zero(factor.rows())).all() ) {
      return desired_resolution;
    } else {
      return desired_resolution + factor - rem;
    }
  }


  FoProperties FoProperties::expand(ArrayFi factor) const
  {
    return FoProperties(m_output_domain.expand(factor),
        m_input_domain.expand(factor));
  }

  FoProperties FoProperties::expandOutputTo(ArrayFi cluster_shape) const
  {
    assert( (cluster_shape.binaryExpr(m_output_domain.clusterShape(),
            std::modulus<int>())).isZero() );

    ArrayFi factor = cluster_shape / m_output_domain.clusterShape();
    return expand(factor);
  }

  FoProperties FoProperties::expandInputTo(ArrayFi cluster_shape) const
  {
    assert( (cluster_shape.binaryExpr(m_input_domain.clusterShape(),
            std::modulus<int>())).isZero() );

    ArrayFi factor = cluster_shape / m_input_domain.clusterShape();
    return expand(factor);

  }

}

