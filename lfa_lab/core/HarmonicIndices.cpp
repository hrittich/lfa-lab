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

#include "HarmonicIndices.h"
#include "FrequencyIndices.h"

#include "MathUtil.h"

namespace lfa {

HarmonicIndices::HarmonicIndices(ArrayFi freq_per_dim,
                                 ArrayFi multiplier,
                                 ArrayFi shift)
 :  m_freq_per_dim(freq_per_dim),
    m_multiplier(multiplier),
    m_shift(shift)
{

}

vector<ArrayFi> HarmonicIndices::elements() const
{
    FrequencyIndices base_set(m_freq_per_dim);
    vector<ArrayFi> base_es = base_set.elements();

    vector<ArrayFi> es(base_es.size());
    for (size_t i = 0; i < base_es.size(); ++i) {
        es[i] = m_multiplier * base_es[i] + m_shift;
    }

    return es;
}

vector<int> HarmonicIndices::indices() const
{
    vector<ArrayFi> es = elements();
    vector<int> inds(es.size());

    FrequencyIndices super_set(m_freq_per_dim * m_multiplier);

    for (size_t i = 0; i < es.size(); ++i) {
        inds[i] = super_set.indexOf(es[i]);
    }

    return inds;
}

ArrayFd HarmonicIndices::baseFrequency(ArrayFd base, ArrayFd step_size) const
{
    return base + 2 * pi * m_shift.cast<double>() /
        (step_size * (m_freq_per_dim * m_multiplier).cast<double>());
}

bool HarmonicIndices::operator== (const HarmonicIndices& other) const
{
    return
        (m_freq_per_dim == other.m_freq_per_dim).all() &&
        (m_multiplier == other.m_multiplier).all() &&
        (m_shift == other.m_shift).all();
}

}

