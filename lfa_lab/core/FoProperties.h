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

#ifndef LFA_FO_PROPERTIES_H
#define LFA_FO_PROPERTIES_H

#include "Common.h"
#include "SplitFrequencyDomain.h"

namespace lfa {

  class FoProperties {
    public:
      FoProperties(
          SplitFrequencyDomain output_domain = SplitFrequencyDomain(),
          SplitFrequencyDomain input_domain = SplitFrequencyDomain());

      FoProperties operator+ (const FoProperties& other) const;
      FoProperties operator- (const FoProperties& other) const {
        return (*this) + other;
      }

      friend FoProperties operator* (complex<double> scalar,
                                     const FoProperties& self)
      {
        return self;
      }

      FoProperties operator* (const FoProperties& other) const;

      FoProperties inverse() const;
      FoProperties adjoint() const;

      Grid outputGrid() const { return m_output_domain.grid(); }
      Grid inputGrid() const { return m_input_domain.grid(); }

      /** Compute the smallest fine grid resolution larger than
       * desired_resolution that can sample the operator.
       */
      ArrayFi adjustResolution(ArrayFi desired_resolution) const;

      /** The output domain. */
      SplitFrequencyDomain output() const { return m_output_domain; }
      /** The input domain. */
      SplitFrequencyDomain input() const { return m_input_domain; }

      int dimension() const { return m_output_domain.dimension(); }

      /** Are the diagonal block of the block matrix representation
       * rectangular? */
      bool rectangular() const { return output() == input(); }

      FoProperties expand(ArrayFi factor) const;

      FoProperties expandOutputTo(ArrayFi cluster_shape) const;
      FoProperties expandInputTo(ArrayFi cluster_shape) const;
    private:
      SplitFrequencyDomain m_output_domain;
      SplitFrequencyDomain m_input_domain;
  };

}

#endif
