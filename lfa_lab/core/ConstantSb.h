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

#ifndef LFA_CONSTANT_SB_H
#define LFA_CONSTANT_SB_H

#include "Common.h"
#include "SymbolBuilder.h"
#include "ClusterSymbol.h"

namespace lfa {

  class ConstantSb : public SymbolBuilder
  {
    public:
      ConstantSb();
      ConstantSb(ClusterSymbol symbol, Grid output_grid, Grid input_grid);

      FoProperties properties();
      Symbol generate(const SamplingProperties& conf);
    private:
      ClusterSymbol m_symbol;
      Grid m_output_grid;
      Grid m_input_grid;
  };

  ClusterSymbol flat_restriction_cluster_symbol(Grid output_grid, Grid input_grid);
  ConstantSb flat_restriction_sb(Grid output_grid, Grid input_grid);

  ClusterSymbol flat_interpolation_cluster_symbol(Grid output_grid, Grid input_grid);
  ConstantSb flat_interpolation_sb(Grid output_grid, Grid input_grid);
}

#endif
