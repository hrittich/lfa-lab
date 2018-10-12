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

  vim: set filetype=cpp:
*/

%module extension

%{
#define SWIG_FILE_WITH_INIT
#include <lfa_lab/core/lfa.h>
#include <sstream>

using namespace lfa;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
%}

%include "std_string.i"
%include "std_complex.i"
// %include "exception.i"

%init %{
    // load NumPy (IMPORTANT!!!)
    import_array();
%}

// enable exception handling
%include exception.i
%exception {
    try {
        $action
    } catch(const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

%include "complex.i"
%include "Eigen.i"
%include "Symbol.i"
%include "shared_ptr.i"
%include "FoStencil.i"
%include "FoProperties.i"
%include "SamplingProperties.i"
%include "SplitFrequencyDomain.i"
%include "SystemSymbol.i"
%include "BlockSb.i"
%include "NdArray.i"
%include "NdRange.i"
%include "SparseStencil.i"
%include "DenseStencil.i"
%include "Grid.i"
%include "ConstantSb.i"
%include "BlockStencil.i"
%include "HpFilterSb.i"
%include "SystemSymbolProperties.i"
%include "BdMatrix.i"

// =========================================================

DenseStencil stencil_poisson2d(ArrayFd h, double eps = 1.0);


//%rename("%(strip:[mk])s") "";

//SfoRef mkFoColJac(const FoSplittingSystem& L, double weight = 1.0);

BlockStencil diff_jump_stride_2d(ArrayFd h,
                                 int blocksize,
                                 double density1,
                                 double density2);
DenseStencil ml_interpolation_stencil(int dim);
DenseStencil fw_restriction(int dim);
BlockStencil flux_conserving_int_2d(BlockStencil L);

void enable_fpe();

