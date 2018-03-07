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

#include "FrequencyIndices.h"
#include "MathUtil.h"
#include "HarmonicClusters.h"
#include "HarmonicIndices.h"

#include <iterator>
#include <algorithm>
#include <stdexcept>

namespace lfa {

using std::back_inserter;

FrequencyIndices::FrequencyIndices(ArrayFi freq_per_dim)
 :  m_freq_per_dim(freq_per_dim)
{
}


vector<HarmonicIndices> FrequencyIndices::split(ArrayFi nsubsets)
{
    assert((mod(m_freq_per_dim, nsubsets) == ArrayFi::Zero(dimension())).all());

    ArrayFi subset_freq_per_dim = m_freq_per_dim / nsubsets;

    vector<HarmonicIndices> subsets;
    NdRange shifts(nsubsets);
    for (NdRange::iterator it = shifts.begin(); it != shifts.end(); ++it) {
        subsets.push_back(HarmonicIndices(subset_freq_per_dim, nsubsets, *it));
    }

    return subsets;
}

vector<HarmonicIndices> FrequencyIndices::splitIntoSizeOf(ArrayFi subset_freq_per_dim)
{
    assert((mod(m_freq_per_dim, subset_freq_per_dim) == ArrayFi::Zero(dimension())).all());

    ArrayFi nsubsets = m_freq_per_dim / subset_freq_per_dim;
    return split(nsubsets);
}

vector<ArrayFi> FrequencyIndices::elements() const
{
    NdRange g(m_freq_per_dim);

    vector<ArrayFi> es;
    copy(g.begin(), g.end(), back_inserter(es));

    return es;
}

vector<ArrayFd> FrequencyIndices::frequencies(Frequency base, ArrayFd step_size) const
{
    // TODO

    assert( (base >= Frequency::Zero(dimension())).all() );
    assert( (base < 2 * pi / step_size * m_freq_per_dim.cast<double>()).all() );

    vector<ArrayFd> shs = shifts(step_size);
    vector<ArrayFd> freqs(shs.size());

    for (size_t i = 0; i < shs.size(); ++i) {
        freqs[i] = base + shs[i];
    }

    return freqs;
}

vector<ArrayFd> FrequencyIndices::shifts(ArrayFd step_size) const
{
    vector<ArrayFi> es = elements();
    vector<ArrayFd> shs(es.size());

    for (size_t i = 0; i < es.size(); ++i) {
        shs[i] = 2 * pi * es[i].cast<double>() /
                    (step_size * m_freq_per_dim.cast<double>());
    }

    return shs;
}

}
