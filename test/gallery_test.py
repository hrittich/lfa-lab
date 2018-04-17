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


from unittest import TestCase

import lfa_lab
from lfa_lab.gallery import *


class GalleryTest(TestCase):

    def test_ml_interpolation_2(self):
        grid = lfa_lab.Grid(2)
        s = ml_interpolation_stencil(grid)

        expect = [ 1/4, 2/4, 1/4,
                   2/4, 4/4, 2/4,
                   1/4, 2/4, 1/4 ]

        for (o, v), v_expect in zip(s, expect):
            self.assertAlmostEqual(v, v_expect)

    def test_ml_interpolation_3(self):
        grid = lfa_lab.Grid(3)
        s = ml_interpolation_stencil(grid)

        expect = [
            1.0/8, 1.0/4, 1.0/8,
            1.0/4, 1.0/2, 1.0/4,
            1.0/8, 1.0/4, 1.0/8,

            1.0/4, 1.0/2, 1.0/4,
            1.0/2, 1.0  , 1.0/2,
            1.0/4, 1.0/2, 1.0/4,

            1.0/8, 1.0/4, 1.0/8,
            1.0/4, 1.0/2, 1.0/4,
            1.0/8, 1.0/4, 1.0/8 ]


        for (o, v), v_expect in zip(s, expect):
            self.assertAlmostEqual(v, v_expect)


    def test_fw_restriction_2(self):
        grid = lfa_lab.Grid(2)
        s = fw_restriction_stencil(grid)

        expect = [
            0.0625, 0.125, 0.0625,
            0.125,  0.25,  0.125,
            0.0625, 0.125, 0.0625 ]

        for (o, v), v_expect in zip(s, expect):
            self.assertAlmostEqual(v, v_expect)

    def test_fw_restriction_3(self):
        grid = lfa_lab.Grid(3)
        s = fw_restriction_stencil(grid)

        expect = [
            1.0/64, 2.0/64, 1.0/64,
            2.0/64, 4.0/64, 2.0/64,
            1.0/64, 2.0/64, 1.0/64,

            2.0/64, 4.0/64, 2.0/64,
            4.0/64, 8.0/64, 4.0/64,
            2.0/64, 4.0/64, 2.0/64,

            1.0/64, 2.0/64, 1.0/64,
            2.0/64, 4.0/64, 2.0/64,
            1.0/64, 2.0/64, 1.0/64 ]

        for (o, v), v_expect in zip(s, expect):
            self.assertAlmostEqual(v, v_expect)





