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

#ifndef LFA_STENCIL_GALLERY_H
#define LFA_STENCIL_GALLERY_H

#include "Common.h"
#include <map>
#include "SystemStencil.h"
#include "BlockStencil.h"

namespace lfa {

  DenseStencil stencil_poisson2d(ArrayFd h, double eps = 1.0);
  DenseStencil stencil_poisson3d(ArrayFd h, ArrayFd eps = ArrayFd::Ones(3));

  SystemStencil stencil_biharmonic_2d(ArrayFd h, ArrayFd);

  DenseStencil stencil_smooth(int distance);
  complex<double> symbol_smooth(int distance, ArrayFd t, ArrayFd h);

  BlockStencil diff_jump_stride_2d(Array2d h,
                                   int blocksize,
                                   double density1,
                                   double density2);

  BlockStencil flux_conserving_int_2d(BlockStencil L);

  DenseStencil ml_interpolation_stencil(int dim);
  DenseStencil fw_restriction(int dim);
}

#endif
