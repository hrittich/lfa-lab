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

#include <gtest/gtest.h>
#include "DenseStencil.h"
#include "Stencil2d.h"
#include "BlockStencil.h"

using namespace lfa;

TEST(Stencil, Consistency)
{
    Vector2i start(-1, -1);
    Vector2i end(1, 1);

    DenseStencil L(start, end);

    EXPECT_EQ( 2, L.dimension() );
}

TEST(DenseStencil, SimpleStoreAndRead)
{
    Vector2i start(-1, -1);
    Vector2i end(1, 1);

    DenseStencil L(start, end);

    Vector2i pos = Vector2i::Zero();
    for (pos(1) = -1; pos(1) <= 1; ++pos(1))
    {
        for (pos(0) = -1; pos(0) <= 1; ++pos(0))
        {
            L(pos) = (pos(1) + 1) * pos(0);
        }
    }

    for (pos(1) = -1; pos(1) <= 1; ++pos(1))
    {
        for (pos(0) = -1; pos(0) <= 1; ++pos(0))
        {
            EXPECT_EQ( complex<double>((pos(1) + 1) * pos(0)), L(pos) );
        }
    }

    DenseStencil::Iterator it(L);
    for (pos(1) = -1; pos(1) <= 1; ++pos(1))
    {
        for (pos(0) = -1; pos(0) <= 1; ++pos(0))
        {
            EXPECT_EQ( complex<double>((pos(1) + 1) * pos(0)), it.value() );
            EXPECT_EQ( true, (bool)it );
            ++it;
        }
    }
    EXPECT_FALSE( it );
}

TEST(DenseStencil, Multiply)
{
    Vector2i start(-1, -1);
    Vector2i end(1, 1);

    DenseStencil L(start, end);

    L(Vector2i(0, -1)) = -1;
    L(Vector2i(-1, 0 )) = -1;
    L(Vector2i(0,0)) = 4;
    L(Vector2i(1, 0)) = -1;
    L(Vector2i(0, 1)) = -1;

    DenseStencil R = L*L;

    EXPECT_EQ(complex<double>(20), R(Vector2i(0,0)));
    EXPECT_EQ(complex<double>(-8), R(Vector2i(0,1)));
    EXPECT_EQ(complex<double>(2), R(Vector2i(1,1)));
    EXPECT_EQ(complex<double>(1), R(Vector2i(0,2)));
    EXPECT_EQ(complex<double>(0), R(Vector2i(2,2)));
}

TEST(DenseStencil, center)
{
    Vector2i start(-1, -1);
    Vector2i end(1, 1);

    DenseStencil L(start, end);

    DenseStencil::Iterator it(L);
    for (int i = 0; i < 4; ++it, ++i)
    {
        EXPECT_TRUE(it.beforeCenter());
        EXPECT_FALSE(it.atCenter());
    }

    EXPECT_FALSE(it.beforeCenter());
    EXPECT_TRUE(it.atCenter());
    ++it;

    for (int i = 0; i < 4; ++it, ++i)
    {
        EXPECT_FALSE(it.beforeCenter());
        EXPECT_FALSE(it.atCenter());
    }
}

TEST(DenseStencil, lower)
{
    Stencil2d L(-1, -1, 1, 1);

    L(-1, -1) = 4;
    L(0, -1) = 4;
    L(1, -1) = 4;
    L(-1, 0) = 4;
    L(0, 0) = 1;
    L(1, 0) = 1;
    L(-1, 1) = 1;
    L(0, 1) = 1;
    L(1, 1) = 1;

    int count = 0;

    DenseStencil lower = L.lower();
    for (DenseStencil::Iterator it(lower); it; ++it) {
        EXPECT_TRUE(it.value() == complex<double>(0)
                    || it.value() == complex<double>(4));
        if (it.value() == complex<double>(4))
            count++;
    }

    EXPECT_EQ(count, 4);
}

TEST(DenseStencil, BlockStencil)
{
    Vector2i asize(2, 2);

    BlockStencil L(asize);

    Stencil2d l(-1, -1, 1, 1);

    l(0, -1) = -1;
    l(-1, 0) = -1;
    l(0,0 ) = 4;
    l(1, 0) = -1;
    l(0, 1) = -1;

    L(Vector2i(1,0)) = l;
    L(Vector2i(0,1)) = l;

    EXPECT_EQ(L.frequencyRange(), 1);


    BlockStencil D = L;
}

TEST(DenseStencil, BlockStencilMultiply)
{
    BlockStencil I( Array2i(2, 2) );
    BlockStencil A( Array2i(2, 2) );

    Stencil2d eye(0, 0, 0, 0);
    eye(0,0) = 1;

    for (BlockStencil::Iterator it(I); it; ++it) {
        I(it.pos()) = eye;
    }

    for (BlockStencil::Iterator it(A); it; ++it) {

        DenseStencil a( Array2i(1, 1), Array2i(3, 3) );

        for (DenseStencil::Iterator inner_it(a); inner_it; ++inner_it) {
            a(inner_it.pos()) = rand();
        }

        A(it.pos()) = a;
    }

    // check that the multiplication is invariant under the identity
    BlockStencil C = A * I;

    EXPECT_TRUE( (C.shape() == A.shape()).all() );

    for (BlockStencil::Iterator it(A); it; ++it) {
        const DenseStencil& c = C(it.pos());
        const DenseStencil& a = A(it.pos());

        EXPECT_TRUE( (c.startIndex() == a.startIndex()).all() );
        EXPECT_TRUE( (c.endIndex() == a.endIndex()).all() );

        for (DenseStencil::ConstIterator inner_it(a); inner_it; ++inner_it) {
            EXPECT_EQ(c(inner_it.pos()), a(inner_it.pos()));
        }
    }
    EXPECT_EQ(C, A);

    // the other way round
    C = I * A;

    for (BlockStencil::Iterator it(A); it; ++it) {
        const DenseStencil& c = C(it.pos());
        const DenseStencil& a = A(it.pos());

        EXPECT_TRUE( (c.startIndex() == a.startIndex()).all() );
        EXPECT_TRUE( (c.endIndex() == a.endIndex()).all() );

        for (DenseStencil::ConstIterator inner_it(a); inner_it; ++inner_it) {
            EXPECT_EQ(c(inner_it.pos()), a(inner_it.pos()));
        }
    }
    EXPECT_EQ(C, A);


    // individual 1D test
    {
        ArrayFi v = ArrayFi::Ones(1);
        DenseStencil a1( -1*v, 0*v );
        a1(-1*v) = 2;
        a1( 0*v) = 1;

        DenseStencil a2( ArrayFi::Zero(1), ArrayFi::Zero(1) );
        a2(0*v) = 1;

        DenseStencil b1( -1*v, 1*v);
        b1(-1*v) = 1;
        b1( 0*v) = 1;
        b1( 1*v) = 1;

        DenseStencil b2( -1*v, 1*v);
        b2( 0*v) = -1;
        b2(-1*v) = -1;
        b2( 1*v) = -1;

        BlockStencil A(2*v), B(2*v);
        A(0*v) = a1;
        A(1*v) = a2;

        B(0*v) = b1;
        B(1*v) = b2;

        BlockStencil C = A*B;

        EXPECT_EQ( C(0*v)(-2*v), complex<double>(-2));
        EXPECT_EQ( C(0*v)(-1*v), complex<double>(-1));
        EXPECT_EQ( C(0*v)(-0*v), complex<double>(-1));
        EXPECT_EQ( C(0*v)( 1*v),  complex<double>(1));

        EXPECT_EQ( C(1*v)(-1*v), complex<double>(-1));
        EXPECT_EQ( C(1*v)(-0*v), complex<double>(-1));
        EXPECT_EQ( C(1*v)( 1*v), complex<double>(-1));
    }

}

TEST(DenseStencil, BlockStencilAdjoint)
{
    DenseStencil a(Array2i(1,0),Array2i(1,0)),
            b(Array2i(0,1),Array2i(0,1));

    a(Array2i(1,0)) = 1;
    b(Array2i(0,1)) = 2;

    BlockStencil L(Array2i(2,2));

    L(Array2i(0,0)) = a;
    L(Array2i(1,0)) = b;
    L(Array2i(0,1)) = b;
    L(Array2i(1,1)) = a;

    BlockStencil T = L.adjoint();

    EXPECT_EQ(T(Array2i(0,0))(Array2i(0, -1)), complex<double>(2));
    EXPECT_EQ(T(Array2i(1,0))(Array2i(-1, 0)), complex<double>(1));
    EXPECT_EQ(T(Array2i(0,1))(Array2i(-1, 0)), complex<double>(1));
    EXPECT_EQ(T(Array2i(1,1))(Array2i(0, -1)), complex<double>(2));
}

