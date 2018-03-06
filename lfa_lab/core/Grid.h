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

#ifndef LFA_FO_GRID_H
#define LFA_FO_GRID_H

#include "Common.h"
#include "FoContext.h"

namespace lfa {

  class Grid {
    public:
      explicit Grid(int dimension = 0,
                    ArrayFd step_size = ArrayFd::Zero(0));
      explicit Grid(shared_ptr<FoContext> ctx);
      explicit Grid(ArrayFi spacing, ArrayFd step_size = ArrayFd::Zero(0));
      explicit Grid(ArrayFi spacing, shared_ptr<FoContext> ctx);

      Grid coarse(ArrayFi factor);

      /** The step size of the grid. */
      ArrayFd step_size() const;

      const ArrayFi& spacing() const { return m_spacing; }

      /** The coarsing factor to the other grid.
       *
       * Warning: Other has to be the coarser grid, i.e.,
       *   other.m_spacing >= this->m_spacing
       */
      ArrayFi coarsening_factor(const Grid& other);

      int dimension() const { return m_spacing.rows(); }

      bool operator== (const Grid& other) const {
        return (m_spacing == other.m_spacing).all();
      }

      ArrayFd finestStepSize() const;

      shared_ptr<FoContext> context() { return m_ctx; }
      shared_ptr<const FoContext> context() const { return m_ctx; }
    private:
      ArrayFi m_spacing;
      shared_ptr<FoContext> m_ctx;
  };

  ostream& operator<< (ostream& os, const Grid& grid);
}

#endif
