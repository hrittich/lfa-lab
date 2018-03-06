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

#include "ExEigenSolver.h"

#include "EigenSolver.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include "Common.h"

using std::strcpy;
using std::stringstream;
using std::endl;
using std::runtime_error;
using std::cerr;
using std::cout;
using std::complex;

extern "C" void znaupd_(
        int* ido,
        const char* bmat,
        int* n,
        const char* which,
        int* nev,
        double* tol,
        complex<double>* resid,
        int* ncv,
        complex<double>* v,
        const int* ldv,
        int* iparam,
        int* ipntr,
        complex<double>* workd,
        complex<double>* workl,
        int* lworkl,
        double* rwork,
        int* info );

extern "C" void zneupd_(
        int* rvec,
        const char* howmny,
        int* select,
        complex<double>* d,
        complex<double>* z,
        const int* ldz,
        complex<double>* sigma,
        complex<double>* workev,
        const char* bmat,
        int* n,
        const char* which,
        int* nev,
        double* tol,
        complex<double>* resid,
        int* ncv,
        complex<double>* v,
        const int* ldv,
        int* iparam,
        int* ipntr,
        complex<double>* workd,
        complex<double>* workl,
        int* lworkl,
        double* rwork,
        int* info );

namespace lfa {

std::complex<double> eigenvalue_max_magnitude(const MatrixXcd& M)
{
#ifdef WITH_ARPACK
    if (M.rows() <= 32)
    {
#endif
        // "directly" solve the eigenvalue problem, if ARPACK is not used or
        // the problem is small enough
        VectorXcd eigs = eigenvalues(M);

        double max_amp = abs(eigs(0));
        for (int i = 1; i < M.rows(); ++i)
        {
            max_amp = std::max(max_amp, abs(eigs(i)));
        }

        //std::cerr << eigs.transpose().array().abs() << std::endl;

        return max_amp;

#ifdef WITH_ARPACK
    }

    int ido = 0; // first call to RCI
    char bmat = 'I'; // standard eigenvalue problem A*x = lambda*x
    int n = M.rows(); // Dimension of the eigenproblem.
    const char* which = "LM"; // want the NEV eigenvalues of largest modulus
    int nev   = 1;  // number of eigenvalues approximated
    double tol = 0; // zero means machine precision
    VectorXcd resid(n);
    int ncv = 30; // The number of Lanczos basis vectors to use through

    const int ldv = n;
    MatrixXcd v(ldv, ncv);

    // parameters
    VectorXi iparam = VectorXi::Zero(11);
    int ishfts = 1;  // exact shifts with respect to the current Hessenberg matrix H.
    int maxitr = 300; // maximum number of Arnoldi update iterations allowed.
    int mode   = 1;

    iparam(0) = ishfts;
    iparam(2) = maxitr;
    iparam(6) = mode;

    // Pointer to mark the starting locations in the WORKD and WORKL
    VectorXi ipntr(14);

    // Work arrays
    VectorXcd workd(3*n);

    int lworkl  = 3 * (ncv*ncv) + 5*ncv;
    VectorXcd workl(lworkl);

    VectorXd  rwork(ncv);

    int info = 0; // random start vector

    assert(ncv - nev >= 2);
    assert(ncv <= n);

    /*-------------------------------------------*
     * M A I N   L O O P (Reverse communication) *
     *-------------------------------------------*/

    #pragma omp critical
    {
        while (true)
        {
            /*---------------------------------------------*
             * Repeatedly call the routine CNAUPD and take *
             * actions indicated by parameter IDO until    *
             * either convergence is indicated or maxitr   *
             * has been exceeded.                          *
             *---------------------------------------------*/

            znaupd_(
                &ido,
                &bmat,
                &n,
                which,
                &nev,
                &tol,
                resid.data(),
                &ncv,
                v.data(),
                &ldv,
                iparam.data(),
                ipntr.data(),
                workd.data(),
                workl.data(),
                &lworkl,
                rwork.data(),
                &info );

            if (ido == -1 || ido == 1)
            {
                /*-------------------------------------------*
                 * Perform matrix vector multiplication      *
                 *------------------------------------------*/

                // input:   workd.segment( ipntr(0), n )
                // output:  workd.segment( ipntr(1), n )

                workd.segment( ipntr(1)-1, n ) =
                    M * workd.segment( ipntr(0)-1, n );
            } else {
                break;
            }
        }
    } // omp critical

    /*----------------------------------------*
     * Either we have convergence or there is *
     * an error.                              *
     *----------------------------------------*/

    if (info < 0) {
        /*-------------------------*
        * Error message, check the *
        * documentation in CNAUPD  *
        *--------------------------*/

        stringstream s;
        s << "Error with _naupd, info = " << info
            << " Check the documentation of _naupd" << endl;
        throw runtime_error(s.str());
    } else {

        /*-------------------------------------------*
         * No fatal errors occurred.                 *
         * Post-Process using CNEUPD.                *
         *                                           *
         * Computed eigenvalues may be extracted.    *
         *                                           *
         * Eigenvectors may also be computed now if  *
         * desired.  (indicated by rvec = .true.)    *
         *-------------------------------------------*/

        int rvec = false;
        const char* howmny = "A";
        VectorXi select(ncv);
        VectorXcd d(ncv);
        complex<double> sigma; // not referenced
        VectorXcd workev(3*ncv);
        int ierr;

        #pragma omp critical
        zneupd_ (
            &rvec, howmny, select.data(), d.data(), v.data(), &ldv, &sigma,
            workev.data(), &bmat, &n, which, &nev, &tol, resid.data(), &ncv,
            v.data(), &ldv, iparam.data(), ipntr.data(), workd.data(),
            workl.data(), &lworkl, rwork.data(), &ierr);

        /*----------------------------------------------*
         * Eigenvalues are returned in the one          *
         * dimensional array D.  The corresponding      *
         * eigenvectors are returned in the first NCONV *
         * (=IPARAM(5)) columns of the two dimensional  *
         * array V if requested.  Otherwise, an         *
         * orthogonal basis for the invariant subspace  *
         * corresponding to the eigenvalues in D is     *
         * returned in V.                               *
         *----------------------------------------------*/

        if ( ierr != 0)
        {
            /*------------------------------------*
             * Error condition:                   *
             * Check the documentation of CNEUPD. *
             *------------------------------------*/

            stringstream s;
            switch(ierr) {
                case -14:
                    s << "CNAUPD did not find any eigenvalues to sufficient accuracy.";
                    break;
                default:
                    s << "Error with _neupd, info = " << ierr <<
                        " Check the documentation of _neupd.";
            }
            throw runtime_error(s.str());
        }

        /*-------------------------------------------%
         | Print additional convergence information. |
         %-------------------------------------------*/

        if (info == 1)
        {
            cerr << "Warning: Maximum number of iterations reached." << endl;
        } else if (info == 3) {
            cerr << "No shifts could be applied during implicit " <<
                "Arnoldi update, try increasing NCV" << endl;
        }
#if 0
        cout <<
            "_NDRV1" << endl <<
            "====== " << endl <<
            endl <<
            "Size of the matrix is " << n << endl <<
            "The number of Ritz values requested is " << nev << endl <<
            "The number of Arnoldi vectors generated (NCV) is " << ncv << endl <<
            "What portion of the spectrum: " << which << endl <<
            "The number of converged Ritz values is " << nconv << endl <<
            "The number of Implicit Arnoldi update iterations taken is " << iparam[2] << endl <<
            "The number of OP*x is " << iparam[8] << endl <<
            "The convergence criterion is " << tol << endl << endl;

        cout << "Eigenvalues: ";
        for (int i = 0; i < nev; ++i)
        {
            cout << abs(d(i)) << " ";
        }
        cout << endl << endl;
#endif

        return d(nev-1);
    }

#endif

}

}

