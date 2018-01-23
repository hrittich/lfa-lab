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

#include "Grid.h"
#include "MathUtil.h"

namespace lfa {

Grid::Grid(int dimension, ArrayFd step_size)
  : m_spacing(ArrayFi::Ones(dimension)),
    m_ctx(new FoContext(dimension, step_size))
{

}

Grid::Grid(shared_ptr<FoContext> ctx)
  : m_spacing(ArrayFi::Ones(ctx->dimension())),
    m_ctx(ctx)
{

}

Grid::Grid(ArrayFi spacing, ArrayFd step_size)
  : m_spacing(spacing),
    m_ctx(new FoContext(spacing.rows(), step_size))
{

}

Grid::Grid(ArrayFi spacing, shared_ptr<FoContext> ctx)
  : m_spacing(spacing),
    m_ctx(ctx)
{

}

Grid Grid::coarse(ArrayFi factor)
{
    return Grid(factor * m_spacing, context());
}

ArrayFd Grid::step_size() const
{
    return context()->step_size() * m_spacing.cast<double>();
}

ArrayFi Grid::coarsening_factor(const Grid& other)
{
    if ((other.m_spacing.binaryExpr(m_spacing, std::modulus<int>())
            != ArrayFi::Zero(m_spacing.rows())).any()) {
        throw runtime_error("Other grid is not a coarsening of this one.");
    }

    return other.m_spacing / m_spacing;
}

ArrayFd Grid::finestStepSize() const
{
    return context()->step_size();
}

ostream& operator<< (ostream& os, const Grid& grid)
{
    os << "Grid(spacing = " << ListFmt(grid.spacing()) << ")";
    return os;
}


}
