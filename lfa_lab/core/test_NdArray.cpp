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

#include "NdArray.h"

using namespace lfa;

TEST(NdArray, simple)
{
    NdArray<int> a(Array2i(3, 6));

    NdRange g = a.shape();
    for (NdRange::iterator p = g.begin();
            p != g.end(); ++p)
    {
        a(*p) = (*p)(0) + 10 * (*p)(1);
    }

    for (NdRange::iterator p = g.begin();
            p != g.end(); ++p)
    {
        EXPECT_EQ((*p)(0) + 10 * (*p)(1), a(*p));
    }


}
