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

#ifndef LFA_FO_STENCIL2_H
#define LFA_FO_STENCIL2_H

#include "Common.h"
#include "SparseStencil.h"
#include "Grid.h"
#include "Symbol.h"
#include "SamplingProperties.h"
#include "DiscreteDomain.h"
#include "SymbolBuilder.h"

namespace lfa {

class FoStencil2 : public SymbolBuilder {
    public:
        FoStencil2(const SparseStencil& stencil, Grid grid);

        FoProperties properties();

        Symbol generate(const SamplingProperties& conf);
        void fill(SymbolClusterRef cluster, ArrayFi base_index, const DiscreteDomain& domain);
        complex<double> symbolAt(VectorFd frequency);

        int dimension() { return m_grid.dimension(); }
    private:
        SparseStencil m_stencil;
        Grid m_grid;
};

}

#endif
