# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018  Hannah Rittich
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from lfa_lab.util import *

import unittest
from unittest import TestCase

class InterpolationTest(TestCase):

    def test_linear(self):

        self.assertAlmostEqual(
                lin_int((0,0), (1, 1), 0.5), 0.5)
        self.assertAlmostEqual(
                lin_int((0,1), (1, 2), 0.5), 1.5)
        self.assertAlmostEqual(
                lin_int((0,0), (1, 2), 0.25), 0.5)
        self.assertAlmostEqual(
                lin_int((1, 2), (0, 0), 0.25), 0.5)

class ShiftedGridTest(TestCase):

    def test_corners(self):

        grid = ShiftedGrid.from_corners( (-1,), (2, ) )
        expect = [ -1, 0, 1, 2 ]

        for p, e in zip(grid, expect):
            self.assertEqual(p, e)


if __name__ == '__main__':
    unittest.main()

