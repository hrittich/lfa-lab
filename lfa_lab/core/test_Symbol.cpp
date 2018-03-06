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

#include "Common.h"

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include "DenseStencil.h"
#include "FoStencil.h"
#include "Symbol.h"
#include "MathUtil.h"
#include "StencilGallery.h"
#include "BlockSb.h"

using namespace lfa;


TEST(Symbol, Iterator)
{
    HarmonicClusters row_clusters(Array2i(2,1), Array2i(1,2));
    HarmonicClusters col_clusters(Array2i(2,1), Array2i(3,1));

    Symbol sym(row_clusters, col_clusters);

    for (Symbol::iterator iter = sym.begin(); iter != sym.end(); ++iter)
    {
        //cout << iter.base() << endl
        //     << iter.row() << endl
        //     << iter.col() << endl << endl;

        *iter = 1;
    }

    MatrixXcd M(4, 6);
    M << 1, 1, 1, 0, 0, 0,
         1, 1, 1, 0, 0, 0,
         0, 0, 0, 1, 1, 1,
         0, 0, 0, 1, 1, 1;

    EXPECT_LE( (sym.full() - M).norm(), 1e-12 );
}

class Poisson1dFixture : public testing::Test
{
    public:
        Poisson1dFixture();
        double theta(ArrayFi pos);

        Grid fine_grid;
        ArrayFi ones;
        ArrayFi zeros;
        DenseStencil poisson_stencil;
        DenseStencil transport_stencil;
        ArrayFd h;
        ArrayFd base_freq;
        ArrayFi shape;
};

Poisson1dFixture::Poisson1dFixture()
  : fine_grid(1)
{
    ones = ArrayFi::Ones(1);
    zeros = ArrayFi::Zero(1);

    // poisson 1d FD stencil
    poisson_stencil.resize(-ones, ones);
    poisson_stencil(-ones) = -1;
    poisson_stencil(zeros) = 2;
    poisson_stencil(ones)  = -1;

    // transport equation stencil
    transport_stencil.resize(-ones, zeros);
    transport_stencil(-ones) = -1;
    transport_stencil(zeros) = 1;

    h = fine_grid.step_size();

    // the sampling density
    shape.resize(1); shape << 8;

    base_freq = pi / (h * shape.cast<double>());

    shared_ptr<FoContext> ctx(new FoContext(1, h));
    fine_grid = Grid(ctx);
}

double Poisson1dFixture::theta(ArrayFi pos)
{
    return base_freq(0) + pos(0) * (2.0 * pi / (h(0) * shape(0)));
}

class Fixture2d : public testing::Test
{
    public:
        Fixture2d();

        ArrayFi ones, zero;
        ArrayFd h, base_freq;
        ArrayFi shape;
        Grid fine_grid;

        ArrayFd theta(ArrayFi b);
};

Fixture2d::Fixture2d()
  : fine_grid(2)
{
    ones = ArrayFi::Ones(2);
    zero = ArrayFi::Zero(2);

    h = fine_grid.step_size();

    shape.resize(2); shape << 4, 4;

    base_freq = pi / (h * shape.cast<double>());

    shared_ptr<FoContext> ctx(new FoContext(2, h));
    fine_grid = Grid(ctx);
}

ArrayFd Fixture2d::theta(ArrayFi b)
{
    return base_freq + b.cast<double>()
         * (2.0 * pi / (h * shape.cast<double>()));
}



