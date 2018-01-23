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

#include "EigenSolver.h"

#include "Config.h"
#include "MathUtil.h"

#include <iostream>
#include <stdexcept>
using namespace std;

#include <cmath>

#ifdef WITH_LAPACK
extern "C" void zgeev_ (
        const char* JOBVL,
        const char* JOBVR,
        const int* N,
        std::complex<double>* A,  // will be overwritten
        const int* LDA,
        std::complex<double>* W,
        std::complex<double>* VL,
        const int* LDVL,
        std::complex<double>* VR,
        const int* LDVR,
        std::complex<double>* WORK,
        const int* LWORK,
        double* RWORK,
        int* INFO);
#else
    #include <Eigen/Eigenvalues>
#endif

#if defined(WITH_LAPACK) && defined(WITH_OUTER_PARALLEL)
#warning LAPACK is not thread safe. It will be used in a serail fashion. \
    This may degenerate performance.
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace lfa {

Eigen::VectorXcd eigenvalues(const Eigen::MatrixXcd& A)
{
    using Eigen::MatrixXcd;
    using Eigen::VectorXcd;
    using Eigen::VectorXd;

    if (A.rows() != A.cols()) {
        throw logic_error("Expecting a square matrix.");
    }
    assert("Matrix is valid" && is_valid(A));

#ifdef WITH_LAPACK

    int N = A.rows();
    MatrixXcd A_aux = A;
    int LDA = N;
    VectorXcd W(N);
    std::complex<double>* const VL = nullptr;
    int LDVL = 1;
    std::complex<double>* const VR = nullptr;
    int LDVR = 1;

    VectorXcd WORK(1);
    int LWORK = -1;
    VectorXd RWORK(2*N);
    int INFO;

    // gfortran + LAPACK not thread safe
    #pragma omp critical
    {
        // Determine the size of the work array
        zgeev_("N", "N",
                &N, A_aux.data(), &LDA, W.data(),
                VL, &LDVL, VR, &LDVR,
                WORK.data(), &LWORK, RWORK.data(), &INFO);

        LWORK = real(WORK(0));
        WORK.resize(LWORK);

        // do the eigenvalue computation
        zgeev_("N", "N",
                &N, A_aux.data(), &LDA, W.data(),
                VL, &LDVL, VR, &LDVR,
                WORK.data(), &LWORK, RWORK.data(), &INFO);
    }

    if (INFO != 0)
        throw runtime_error("Eigenvalue computation failed.");

    return W;
#else
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigs_comp(A, false);

    return eigs_comp.eigenvalues();
#endif
}

}

