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

#ifndef LFA_STENCIL_2D_H
#define LFA_STENCIL_2D_H

#include "Common.h"
#include <iosfwd>
#include "DenseStencil.h"

namespace lfa {

/** Storage for a 2D stencil. */
class Stencil2d : public DenseStencil
{
    public:
        Stencil2d(const DenseStencil& s);
        Stencil2d(int x0 = 0, int y0 = 0, int x1 = 0, int y1 = 0);

        inline
        complex<double>& operator() (int x, int y) {
            return DenseStencil::operator() (Vector2i(x,y));
        }

        inline complex<double> operator() (int x, int y) const {
            return DenseStencil::operator() (Vector2i(x, y));
        }

        /** Sets the stencil to the given size and fill it with zeros. */
        void resize(int x0, int y0, int x1, int y1);

        inline int x0() const { return startIndex()(0); }
        inline int y0() const { return startIndex()(1); }
        inline int x1() const { return endIndex()(0); }
        inline int y1() const { return endIndex()(1); }

        friend std::ostream& operator<< (std::ostream& os, const Stencil2d& s);
};

}

#endif
