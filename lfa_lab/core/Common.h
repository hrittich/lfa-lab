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

#ifndef LFA_COMMON_H
#define LFA_COMMON_H

#ifndef NDEBUG
#define DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#endif


#include <Config.h>

#include <complex>

#include <iostream>
#include <iosfwd>
#include <stdexcept>
#include <Eigen/Core>

#define LfaArraySize(a) ( sizeof(a) / sizeof(*(a)) )

#ifdef GCC_BOUND_CHECKS
    #include <debug/string>
    #include <debug/vector>
#else
    #include <string>
    #include <vector>
#endif


#ifdef HAVE_STD_SHARED_PTR
    #include <memory>
#else
    #include <boost/shared_ptr.hpp>
#endif

#ifndef HAVE_NULLPTR
#define nullptr NULL
#endif

namespace lfa {

  const size_t MAX_DIMENSION = 5;

  using std::stringstream;
  using std::ostream;
  using std::istream;
  using std::complex;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::logic_error;
  using std::runtime_error;
  using std::out_of_range;
  using std::pair;

#ifdef GCC_BOUND_CHECKS
  using __gnu_debug::vector;
  using __gnu_debug::string;
#else
  using std::vector;
  using std::string;
#endif

#ifdef HAVE_STD_SHARED_PTR
  using std::shared_ptr;
#else
  using boost::shared_ptr;
#endif

  using namespace Eigen;

  /* Fixed Matrix sizes. */
  typedef Matrix<int, Dynamic, 1, ColMajor|AutoAlign, MAX_DIMENSION, 1> VectorFi;
  typedef Matrix<double, Dynamic, 1, ColMajor|AutoAlign, MAX_DIMENSION, 1> VectorFd;
  typedef Matrix< std::complex<double>, Dynamic, 1, ColMajor|AutoAlign, MAX_DIMENSION, 1> VectorFcd;
  typedef Array<int, Dynamic, 1, ColMajor|AutoAlign, MAX_DIMENSION, 1> ArrayFi;
  typedef Array<double, Dynamic, 1, ColMajor|AutoAlign, MAX_DIMENSION, 1> ArrayFd;

  class HarmonicIndices;

  typedef ArrayFd Frequency;

  /** Enable floating point exception. */
  void enable_fpe();
}

namespace std {
  extern template class vector<lfa::ArrayFi>;
  extern template class vector<lfa::ArrayFd>;
  extern template class vector<int>;
  extern template class vector<double>;
}

#endif
