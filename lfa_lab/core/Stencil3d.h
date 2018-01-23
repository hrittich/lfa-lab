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

#ifndef LFA_STENCIL_3D_H
#define LFA_STENCIL_3D_H

#include "Common.h"
#include "DenseStencil.h"

namespace lfa {

/** Storage for a 3D stencil. */
class Stencil3d : public DenseStencil
{
    public:
        Stencil3d(int x0, int y0, int z0,
                int x1, int y1, int z1)
            : DenseStencil(Vector3i(x0, y0, z0), Vector3i(x1, y1, z1) )
        { }

        inline double& operator() (int x, int y, int z) {
            return DenseStencil::operator() ( Vector3i(x,y,z) );
        }

        inline double operator() (int x, int y, int z) const {
            return DenseStencil::operator() ( Vector3i(x,y,z) );
        }
};

}

#endif
