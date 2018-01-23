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

#include "Stencil2d.h"

#include <iomanip>
#include <algorithm>
#include "MathUtil.h"

namespace lfa {

Stencil2d::Stencil2d(const DenseStencil& s)
{
    assert(s.dimension() == 2);
    static_cast<DenseStencil&>(*this) = s;
}

Stencil2d::Stencil2d(int x0, int y0, int x1, int y1)
{
    resize(x0, y0, x1, y1);
}

void Stencil2d::resize(int x0, int y0, int x1, int y1)
{
    Vector2i start(x0, y0), end(x1, y1);
    DenseStencil::resize(start, end);
}

std::ostream& operator<< (std::ostream& os, const Stencil2d& s)
{
    using namespace std;

    for (int y = s.y0(); y <= s.y1(); ++y) {
        for (int x = s.x0(); x <= s.x1(); ++x) {
            os << setw(10) << s(x, y);
        }
        os << endl;
    }

    return os;
}

}

