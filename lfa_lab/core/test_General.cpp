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

#include "MathUtil.h"
#include "NdRange.h"

using namespace lfa;

#include <iterator>
#include <algorithm>

TEST(General, log2)
{
    EXPECT_EQ(0, log2_rn(1));
    EXPECT_EQ(1, log2_rn(2));
    EXPECT_EQ(2, log2_rn(4));
    EXPECT_EQ(3, log2_rn(8));

    EXPECT_TRUE(isPow2(1));
    EXPECT_TRUE(isPow2(2));
    EXPECT_TRUE(isPow2(4));
    EXPECT_TRUE(isPow2(8));

    EXPECT_FALSE(isPow2(3));
    EXPECT_FALSE(isPow2(5));
    EXPECT_FALSE(isPow2(6));
    EXPECT_FALSE(isPow2(7));
}

TEST(General, grid)
{
    using namespace std;
    using namespace Eigen;

    int x[] = { 0, 1, 0, 1, 0, 1, 99 };
    int y[] = { 0, 0, 1, 1, 2, 2, 99 };

    NdRange g(Array2i(2, 3));

    // copy(g.begin(), g.end(), ostream_iterator<ArrayFi>(cout, "\n\n"));

    int i = 0;
    for (NdRange::iterator it = g.begin(); it != g.end(); ++it, ++i) {
        EXPECT_EQ( (*it)(0), x[i] );
        EXPECT_EQ( (*it)(1), y[i] );
    }

    EXPECT_EQ(i, 6);
}

