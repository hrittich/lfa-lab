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

from lfa_lab import *
from lfa_lab.smoother import *
from lfa_lab.gallery import *

import unittest

class SmootherRun(unittest.TestCase):
    """Run all smoothers to see that the code runs."""

    def test_run_jacobi(self):
        fine = Grid(2)
        L = poisson_2d(fine)
        S = jacobi(L)

    def test_run_gs_lex(self):
        fine = Grid(2)
        L = poisson_2d(fine)
        S = gs_lex(L)

    def test_run_rb_jacobi(self):
        g = Grid(2)
        s = poisson_2d(g)
        rb_jacobi(s)

if __name__ == '__main__':
    unittest.main()
