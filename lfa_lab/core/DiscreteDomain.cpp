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

#include "DiscreteDomain.h"
#include "Math.h"

namespace lfa {

DiscreteDomain::DiscreteDomain(
                const SplitFrequencyDomain& domain,
                const SamplingProperties& conf)
  : m_domain(domain),
    m_conf(conf)
{


}

HarmonicClusters DiscreteDomain::harmonics() const
{
    return m_domain.harmonics(resolution());
}

ArrayFd DiscreteDomain::frequency(ArrayFi base_index, ArrayFi cluster_index) const
{
    HarmonicClusters cluster = harmonics();

    ArrayFi g = cluster.globalIndex(base_index, cluster_index);
    return frequency(g);
}


ArrayFd DiscreteDomain::frequency(ArrayFi global_index) const
{
    ArrayFd base_freq = m_conf.base_frequency();
#ifdef DEBUG
    ArrayFd max_base = 2.0 * pi / (step_size() * resolution().cast<double>());
#endif
    assert( (base_freq >= ArrayFd::Zero(dimension())).all() );
    assert( (base_freq <= max_base).all() );


    return base_freq + global_index.cast<double>() * 2.0 * pi
            / (step_size() * resolution().cast<double>());
}

ArrayFd DiscreteDomain::step_size() const
{
    return m_domain.step_size();
}

ArrayFi DiscreteDomain::resolution() const
{
    return m_domain.resolution(m_conf.finest_resolution());
}




}

