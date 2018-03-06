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

#ifndef BLOCK_SB_H
#define BLOCK_SB_H

#include "Common.h"
#include "SymbolBuilder.h"
#include "NdArray.h"

namespace lfa {

  /** Block symbol builder. Combines periodic scalar symbols into a matrix
   * symbol.
   */
  class BlockSb : public SymbolBuilder {
    public:
      BlockSb(Grid grid, ArrayFi period);

      /** Set the scalar symbols that should be combined. */
      void scalarSymbols(NdArray<Symbol> scalars);

      /** The properties of the whole operator. */
      FoProperties properties();

      /** The input and output domain of the scalar operators. */
      SplitFrequencyDomain scalar_domain();

      /** The properties of the scalar operator. */
      FoProperties scalar_properties();

      Symbol generate(const SamplingProperties& conf);

      int dimension() { return m_grid.dimension(); }
    private:
      Grid m_grid;
      ArrayFi m_period;
      NdArray<Symbol> m_scalars;
  };


}

#endif
