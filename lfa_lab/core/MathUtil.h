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

#ifndef LFA_MATH_HELPER_H
#define LFA_MATH_HELPER_H

#include <lfa_lab/core/Common.h>
#include <lfa_lab/core/NdArray.h>

/** @addtogroup Math
 * @{
 */

namespace lfa {

  static const double pi = 3.14159265358979323846;

  /** Integer division with rounding towards zero. */
  int div_rz(int a, int b);

  bool is_valid(std::complex<double> a);
  bool is_valid(double a);

  /** Check if a matrix contains only valid entries, i.e. no NaNs and no INFs.
   * */
  template <typename Derived>
  bool is_valid(const Eigen::DenseBase<Derived>& _A)
  {
    typename Eigen::DenseBase<Derived>::EvalReturnType A = _A.eval();

    for (int j = 0; j < A.cols(); ++j)
      for (int i = 0; i < A.rows(); ++i)
      {
        if ( !is_valid(A.derived()(i,j) ) )
        {
          return false;
        }
      }
    return true;
  }

  /** The wave function. */
  std::complex<double> phi(const ArrayFd& theta, const ArrayFi& pos);

  /** Positive sign. 1 for x > 0, -1 otherwise. */
  inline double psign(double x) {
    if (x >= 0) {
      return 1;
    } else {
      return -1;
    }
  }

  /** Base 2 logarithm, round to negative */
  inline size_t log2_rn(size_t x) {
    size_t l = 0;

    while (x >>= 1) {
      ++l;
    }
    return l;
  }

  /** Is a power of two. */
  inline bool isPow2(size_t x) {

    size_t l = log2_rn(x);
    return (x == (1lu << l));
  }

  /** Mathematical modulus. */
  int mod(int p, int q);
  ArrayFi mod(const ArrayFi& p, const ArrayFi& q);

  int cmod(int p, int q);
  ArrayFi cmod(ArrayFi p, ArrayFi q);

  double cfmod(double p, double q);
  ArrayFd cfmod(const ArrayFd& p, const ArrayFd& q);


  ArrayFi gcd(ArrayFi a, ArrayFi b);
  ArrayFi lcm(ArrayFi a, ArrayFi b);


  bool is_similar(double a, double b);
  bool is_similar_c(complex<double> a, complex<double> b);
  bool is_similar_mc(MatrixXcd A, MatrixXcd B);



  double float_hash(const vector<int>& data, int seed = 0);

  class ListFmt {
    public:
      ListFmt(ArrayFi data);

      friend ostream& operator<< (ostream& os, const ListFmt& f);
    private:
      ArrayFi m_data;
  };

  class IsApprox {
    public:
      IsApprox(double precision);

      bool operator() (double a, double b);
      bool operator() (complex<double> a, complex<double> b);

    private:
      double m_precision;
  };

  ArrayFi lcm(ArrayFi a, ArrayFi b);

  /** Square a number. */
  template <typename number_t>
    double sq(number_t x) { return x*x; }

  /** The square of the absolute value of a complex number. */
  double abs_sq(complex<double> x);

  MatrixXcd to_matrix(NdArray<double> A);
  VectorXcd to_vector(NdArray<double> A);


/**
 * @}
 */


}

#endif
