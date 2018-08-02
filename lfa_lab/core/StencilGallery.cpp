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

#include "StencilGallery.h"

#include "Stencil2d.h"
#include "Stencil3d.h"

#include <iostream>
#include <stdexcept>

namespace lfa {

DenseStencil stencil_poisson2d(ArrayFd h, double eps)
{
    Stencil2d L(-1,-1,1,1);
    double h0 = h(0), h1 = h(1);

    L( 0, -1) = -1 / (h1*h1);
    L(-1,  0) = -1 / (h0*h0) * eps;
    L( 0,  0) =  2 / (h0*h0) * eps + 2 / (h1*h1);
    L( 1,  0) = -1 / (h0*h0) * eps;
    L( 0,  1) = -1 / (h1*h1);

    return L;
}

DenseStencil stencil_poisson3d(ArrayFd h, ArrayFd eps)
{
    Stencil3d L(
        -1, -1, -1,
        1,  1,  1);

    L( 1, 0, 0) = -1.0 / (h(0)*h(0))*eps(0);
    L( 0, 1, 0) = -1.0 / (h(1)*h(1))*eps(1);
    L( 0, 0, 1) = -1.0 / (h(2)*h(2))*eps(2);
    L( 0, 0, 0) =  2.0 / (h(0)*h(0))*eps(0)
                  +2.0 / (h(1)*h(1))*eps(1)
                  +2.0 / (h(2)*h(2))*eps(2);
    L(-1, 0, 0) = -1.0 / (h(0)*h(0))*eps(0);
    L( 0,-1, 0) = -1.0 / (h(1)*h(1))*eps(1);
    L( 0, 0,-1) = -1.0 / (h(2)*h(2))*eps(2);

    return L;
}


SystemStencil stencil_biharmonic_2d(ArrayFd h, ArrayFd)
{
    SystemStencil L(2,2); // 2x2 system

    L(0,0) = stencil_poisson2d(h);
    L(0,1).setZero(2);

    L(1,0).setIdentity(2); L(1,0)( Vector2i::Zero() ) = -1.0;
    L(1,1) = stencil_poisson2d(h);

    return L;
}


DenseStencil stencil_smooth(int distance)
{
    assert(distance > 0);

    ArrayFi ones = ArrayFi::Ones(2);

    DenseStencil s(-distance*ones, distance*ones);

    s(Array2i(-distance, 0)) = 0.25;
    s(Array2i( 0,-distance)) = 0.25;
    s(Array2i( distance, 0)) = 0.25;
    s(Array2i( 0, distance)) = 0.25;

    return s;
}

complex<double> symbol_smooth(int distance, ArrayFd t, ArrayFd h)
{
    assert(distance > 0);

    return 0.5 * cos(t(0)*h(0)*distance)
         + 0.5 * cos(t(1)*h(1)*distance);
}




/*************************************************************************/
/* Stencils for the 2d jumping coefficient problem.                      */
/* See Trottenberg, Oosterlee, Schueller -- Multigrid                    */
/*************************************************************************/
DenseStencil diff_jump_2d(Array2d h,
                     double density_w,
                     double density_e,
                     double density_n,
                     double density_s,
                     double density_c)
{
    if (h(0) != h(1))
        throw std::runtime_error("Only equal step sizes supported.");

    double h_inv = 1/(h(0)*h(1));

    double s_w = -2 * (density_w*density_c)/(density_w + density_c);
    double s_e = -2 * (density_e*density_c)/(density_e + density_c);
    double s_n = -2 * (density_n*density_c)/(density_n + density_c);
    double s_s = -2 * (density_s*density_c)/(density_s + density_c);
    double s_c = -(s_w+s_e+s_n+s_s);

    Stencil2d L(-1, -1, 1, 1);
    L(-1, 0) = s_w * h_inv;
    L(1, 0) = s_e * h_inv;
    L(0, -1) = s_n * h_inv;
    L(0, 1) = s_s * h_inv;
    L(0, 0) = s_c * h_inv;

    return L;
}

double diff_jump_density_at(int blocksize, int x, double density1, double density2)
{
    using namespace std;

    int jump1 = blocksize / 4;
    int jump2 = (3*blocksize) / 4;

    //cerr << "jump at " << jump1 << " and " << jump2 << endl;

    // w density
    if (x < jump1 || x >= jump2) {
        return density1;
    } else {
        return density2;
    }
}

BlockStencil diff_jump_stride_2d(Array2d h,
                                 int blocksize,
                                 double density1,
                                 double density2)
{
    if (blocksize <= 0) {
        throw std::runtime_error("block size to small");
    }

    VectorFi shape = blocksize * VectorFi::Ones(2);
    BlockStencil L(shape);

    for (BlockStencil::Iterator it(L); it; ++it)
    {
        int x = it.pos()(0);
        double density_w = diff_jump_density_at(blocksize, x-1, density1, density2);
        double density_c = diff_jump_density_at(blocksize, x,   density1, density2);
        double density_e = diff_jump_density_at(blocksize, x+1, density1, density2);

        L(it.pos()) = diff_jump_2d(h,
                                   density_w,
                                   density_e,
                                   density_c,
                                   density_c,
                                   density_c);
    }

    return L;
}

/*************************************************************************/



BlockStencil flux_conserving_int_2d(BlockStencil L)
{
    if (L.dimension() != 2)
        throw std::runtime_error("Dimension != 2");

    VectorFi shape = L.shape();
    assert(shape(0) % 2 == 0);
    assert(shape(1) % 2 == 0);


    BlockStencil P(shape);

    DenseStencil Z(Array2i(0,0), Array2i(0,0));
    DenseStencil I(Array2i(0,0), Array2i(0,0));
    I(Array2i(0, 0)) = 1.0;

    // initialize the block stencil with zero
    for (int y = 0; y < shape(1); ++y) {
        for (int x = 0; x < shape(0); ++x) {
            P(Array2i(x,y)) = Z;
        }
    }

    // set the interpolation for the coarse grid points
    for (int y = 0; y < shape(1); y+=2) {
        for (int x = 0; x < shape(0); x+=2) {
            P(Array2i(x,y)) = I;
        }
    }

    // between two points in x-direction
    for (int y = 0; y < shape(1); y+=2) {
        for (int x = 1; x < shape(0); x+=2) {
            const DenseStencil& l = L(Array2i(x,y));

            complex<double>
              center = l(Array2i( 0, -1)) +
                       l(Array2i( 0,  0)) +
                       l(Array2i( 0, +1));
            complex<double>
              left = l(Array2i(-1, -1)) +
                     l(Array2i(-1,  0)) +
                     l(Array2i(-1, +1));
            complex<double>
              right = l(Array2i(+1, -1)) +
                      l(Array2i(+1,  0)) +
                      l(Array2i(+1, +1));

            DenseStencil s(Array2i(-1, 0), Array2i(1, 0));
            s(Array2i(-1, 0)) = -(left / center);
            s(Array2i(1, 0)) = -(right / center);

            P(Array2i(x, y)) = s;
        }
    }

    // between two points in y-direction
    for (int y = 1; y < shape(1); y+=2) {
        for (int x = 0; x < shape(0); x+=2) {
            const DenseStencil& l = L(Array2i(x,y));

            complex<double>
              center = l(Array2i(-1, 0)) +
                       l(Array2i( 0, 0)) +
                       l(Array2i(+1, 0));
            complex<double>
              top = l(Array2i(-1, -1)) +
                    l(Array2i( 0, -1)) +
                    l(Array2i(+1, -1));
            complex<double>
              bottom = l(Array2i(-1, +1)) +
                       l(Array2i( 0, +1)) +
                       l(Array2i(+1, +1));

            DenseStencil s(Array2i(0, -1), Array2i(0, 1));
            s(Array2i(0, -1)) = -(top/center);
            s(Array2i(0, 1)) = -(bottom/center);

            P(Array2i(x, y)) = s;
        }
    }

    // in the center of 4 points
    // indirect interpolation
    for (int y = 1; y < shape(1); y+=2) {
        for (int x = 1; x < shape(0); x+=2) {

            const DenseStencil& l = L(Array2i(x,y));
            complex<double> center = l(Array2i(0,0));

            DenseStencil s(Array2i(-1, -1), Array2i(1, 1));

            s(Array2i(-1,-1)) = - l(Array2i(-1,-1)) / center;
            s(Array2i( 0,-1)) = - l(Array2i( 0,-1)) / center;
            s(Array2i( 1,-1)) = - l(Array2i( 1,-1)) / center;
            s(Array2i(-1, 0)) = - l(Array2i(-1, 0)) / center;
            s(Array2i( 1, 0)) = - l(Array2i( 1, 0)) / center;
            s(Array2i(-1, 1)) = - l(Array2i(-1, 1)) / center;
            s(Array2i( 0, 1)) = - l(Array2i( 0, 1)) / center;
            s(Array2i( 1, 1)) = - l(Array2i( 1, 1)) / center;

            // the values N,S,W,E are unknown thus we use the interpolation in
            // those points to compute the values for those nodes

            s = P.multiplyRight(s, Array2i(x,y));

            P(Array2i(x,y)) = s;
        }
    }

    return P;
}


// full weighting restriction
double restrict_fw2d_data[] = {
    0.0625, 0.125, 0.0625,
    0.125,  0.25,  0.125,
    0.0625, 0.125, 0.0625
};

double restrict_fw3d_data[] = {
    1.0/64, 2.0/64, 1.0/64,
    2.0/64, 4.0/64, 2.0/64,
    1.0/64, 2.0/64, 1.0/64,

    2.0/64, 4.0/64, 2.0/64,
    4.0/64, 8.0/64, 4.0/64,
    2.0/64, 4.0/64, 2.0/64,

    1.0/64, 2.0/64, 1.0/64,
    2.0/64, 4.0/64, 2.0/64,
    1.0/64, 2.0/64, 1.0/64
};


DenseStencil fw_restriction(int dimension)
{
    VectorFi start = VectorFi::Ones(dimension) * (-1);
    VectorFi end = VectorFi::Ones(dimension);

    DenseStencil R(start, end);
    switch(dimension) {
        case 2:
            R.assign(restrict_fw2d_data,
                    restrict_fw2d_data + LfaArraySize(restrict_fw2d_data));
            break;
        case 3:
            R.assign(restrict_fw3d_data,
                    restrict_fw3d_data + LfaArraySize(restrict_fw3d_data));

            break;
        default:
            throw std::runtime_error("unsupported dimension");
    }

    return R;
}

// bilinear interpolation
double interpolate_ml2d_data[] = {
    1.0/4, 2.0/4, 1.0/4,
    2.0/4, 4.0/4, 2.0/4,
    1.0/4, 2.0/4, 1.0/4
};

double interpolate_ml3d_data[] = {
    1.0/8, 1.0/4, 1.0/8,
    1.0/4, 1.0/2, 1.0/4,
    1.0/8, 1.0/4, 1.0/8,

    1.0/4, 1.0/2, 1.0/4,
    1.0/2, 1.0  , 1.0/2,
    1.0/4, 1.0/2, 1.0/4,

    1.0/8, 1.0/4, 1.0/8,
    1.0/4, 1.0/2, 1.0/4,
    1.0/8, 1.0/4, 1.0/8
};

DenseStencil ml_interpolation_stencil(int dimension)
{
    VectorFi start = VectorFi::Ones(dimension) * (-1);
    VectorFi end = VectorFi::Ones(dimension);

    DenseStencil P(start, end);
    switch(dimension) {
        case 2:
            P.assign(interpolate_ml2d_data,
                    interpolate_ml2d_data + LfaArraySize(interpolate_ml2d_data));
            break;
        case 3:
            P.assign(interpolate_ml3d_data,
                    interpolate_ml3d_data + LfaArraySize(interpolate_ml3d_data));

            break;
        default:
            throw std::runtime_error("unsupported dimension");
    }

    return P;
}



}

