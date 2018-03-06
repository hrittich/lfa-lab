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

#include "FoProperties.h"

using namespace lfa;

TEST(FoProperties, Simple)
{
    Grid grid(2); // 1d grid

    SplitFrequencyDomain dom_a(grid, Array2i(3, 3));
    SplitFrequencyDomain dom_b(grid, Array2i(2, 2));
    SplitFrequencyDomain dom_c(grid, Array2i(1, 1));

    FoProperties first(dom_c, dom_a);
    FoProperties second(dom_b, dom_c);

    FoProperties res = first * second;

    EXPECT_EQ(2, res.output().clusterShape()(0));
    EXPECT_EQ(3, res.input().clusterShape()(0));
}
