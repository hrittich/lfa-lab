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

#include "CartIterator.h"

using namespace lfa;

TEST(CartIterator, Simple)
{
    int nums[] = { 0, 1, 2, 3, 4, 5 };

    CartIterator<int*, int*> iter(nums, nums+4, nums, nums+6);
    CartIterator<int*, int*> end(nums, nums+4, nums, nums+6, true);

    int n = 0;
    while (iter != end) {
        // cout << *(iter->first) << " " << *(iter->second) << endl;

        int i = *(iter->first);
        int j = *(iter->second);

        EXPECT_EQ(i + 4*j, n);

        ++iter;
        ++n;
    }

}

