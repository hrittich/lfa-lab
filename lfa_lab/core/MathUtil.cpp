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

#include "MathUtil.h"

#include <cstdlib>

using namespace std;

namespace lfa {

  int div_rz(int a, int b)
  {
    int ma = abs(a);
    int mb = abs(b);

    if ( (a > 0 && b > 0) || (a < 0 && b < 0)) {
      return ma / mb;
    } else {
      return - (ma / mb);
    }
  }

  bool is_valid(double a) {
    return !std::isnan(a) && !std::isinf(a);
  }

  bool is_valid(std::complex<double> a) {
    return is_valid(real(a)) && is_valid(imag(a));
  }

  /** The wave function. */
  std::complex<double> phi(const ArrayFd& theta, const ArrayFi& pos)
  {
    using namespace std;
    using namespace Eigen;

    return exp( complex<double> (0, ArrayFd(theta * pos.cast<double>() ).sum() ) );
  }


  int mod(int p, int q) {
    int r = p % q;
    if (r < 0) {
      return r + q;
    } else {
      return r;
    }
  }

  ArrayFi mod(const ArrayFi& p, const ArrayFi& q) {
    ArrayFi r(p.size());

    for (int i = 0; i < p.size(); ++i) {
      r(i) = mod(p(i), q(i));
    }
    return r;
  }

  int cmod(int p, int q) {
    int r = mod(p, q);
    if (r > q/2) {
      r = r - q;
    }

    return r;
  }

  ArrayFi cmod(ArrayFi p, ArrayFi q)
  {
    ArrayFi r(p.rows());

    for (int i = 0; i < p.size(); ++i) {
      r(i) = cmod(p(i), q(i));
    }
    return r;
  }

  double cfmod(double p, double q)
  {
    double r = fmod(p, q);
    if (r < 0) {
      r += q;
    }

    if (r > q/2) {
      r -= q;
    }
    return r;
  }

  ArrayFd cfmod(const ArrayFd& p, const ArrayFd& q) {
    ArrayFd r(p.size());

    for (int i = 0; i < p.size(); ++i) {
      r(i) = cfmod(p(i), q(i));
    }
    return r;
  }


  int gcd(int a, int b)
  {
    while (b != 0)
    {
      int r = a % b;
      a = b;
      b = r;
    }
    return a;
  }

  ArrayFi gcd(ArrayFi a, ArrayFi b)
  {
    ArrayFi g(a.rows());

    for (int i = 0; i < g.rows(); ++i) {
      g(i) = gcd(a(i), b(i));
    }
    return g;
  }

  template <typename T>
  T lcm_imp(T a, T b)
  {
    T g = gcd(a, b);
    return (a / g) * b;
  }

  ArrayFi lcm(ArrayFi a, ArrayFi b)
  {
    return lcm_imp<ArrayFi>(a, b);
  }


  bool is_similar(double a, double b)
  {
    double scale = std::max(fabs(a), fabs(b));
    scale = max(scale, 1e-292);

    double diff = a-b;
    double adiff = fabs(diff);
    double sdiff = adiff / scale;

    return sdiff < 1e-12;
  }

  bool is_similar_c(complex<double> a, complex<double> b)
  {
    double scale_real = std::max(fabs(real(a)), fabs(real(b)));
    double scale_imag = std::max(fabs(imag(a)), fabs(imag(b)));
    double scale = std::max(scale_real, scale_imag);

    return norm(a - b) / scale < 1e-12;
  }

  bool is_similar_mc(MatrixXcd A, MatrixXcd B)
  {
    double scale = std::max(A.cwiseAbs().maxCoeff(),
        B.cwiseAbs().maxCoeff());

    double diff = (A-B).cwiseAbs().maxCoeff();
    double sdiff = diff / scale;

    return sdiff < 1e-12;
  }


  double float_hash(const vector<int>& data, int seed)
  {
#warning "Will produce different results on 32 and 64 bit machines"
    const int MAX_HASH = 1 << 30;

    int a = 48271;
    int b = seed;
    int c = 12345;
    int sum = 0;

    for (size_t i = 0; i < data.size(); ++i) {
      b = (b + c) * a;
      sum += b * data[i];
    }
    sum = sum % MAX_HASH;

    return static_cast<double>(sum) / (MAX_HASH - 1);
  }

  ListFmt::ListFmt(ArrayFi data)
    : m_data(data)
  { }

  ostream& operator<< (ostream& os, const ListFmt& f)
  {
    os << "[";
    int i = 0;
    while (true) {
      os << f.m_data(i);
      i += 1;

      if (i >= f.m_data.rows())
        break;

      os << ", ";
    }
    os << "]";
    return os;
  }

  IsApprox::IsApprox(double precision)
    : m_precision(precision)
  { }

  bool IsApprox::operator() (double a, double b)
  {
    return fabs(a - b) < m_precision;
  }

  bool IsApprox::operator() (complex<double> a, complex<double> b)
  {
    return abs(a - b) < m_precision;
  }

  double abs_sq(complex<double> x)
  {
    return sq(real(x)) + sq(imag(x));
  }

  MatrixXcd to_matrix(NdArray<double> A) {
    if (A.dimension() != 2) {
      throw logic_error("Only two-dimensional arrays allowed.");
    }

    MatrixXcd result(A.shape()(0), A.shape()(1));
    for (NdRange::iterator p = A.indices().begin();
        p != A.indices().end(); ++p)
    {
      result((*p)(0), (*p)(1)) = A(*p);
    }

    return result;
  }

  VectorXcd to_vector(NdArray<double> A)
  {
    if (A.dimension() != 1) {
      throw logic_error("Only one-dimensional arrays allowed.");
    }

    VectorXcd result(A.shape()(0));
    for (NdRange::iterator p = A.indices().begin();
        p != A.indices().end(); ++p)
    {
      result((*p)(0)) = A(*p);
    }

    return result;
  }

}

