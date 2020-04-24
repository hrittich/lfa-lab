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

import lfa_lab
from lfa_lab import *
import unittest

class OperatorTest(unittest.TestCase):

    def setUp(self):
        self.fine = lfa_lab.Grid(2)
        self.coarse = self.fine.coarse((2,2))

    def test_identity(self):
        I = lfa_lab.operator.identity(self.fine)
        self.assertEqual('id', repr(I))

    def test_zero(self):
        Z = lfa_lab.operator.zero(self.fine)
        self.assertEqual('0', repr(Z))

    def test_sum(self):
        I = lfa_lab.operator.identity(self.fine)
        s = I + I
        self.assertEqual('(+\n  id\n  id)', repr(s))

    def test_product(self):
        I = lfa_lab.operator.identity(self.fine)
        s = I * I
        self.assertEqual('(*\n  id\n  id)', repr(s))

    def test_scalar_product(self):
        I = lfa_lab.operator.identity(self.fine)
        s = 1 * I
        self.assertEqual('(*\n  1\n  id)', repr(s)) 

    def test_stencil(self):
        A = lfa_lab.operator.from_stencil([((0,), 1)], self.fine)
        self.assertEqual('(stencil [((0,), 1.0)])', repr(A))

    def test_hp_filter(self):
        F = operator.hp_filter(self.fine, self.coarse)
        self.assertEqual('(hp_filter\n  (grid 1 1)\n  (grid 2 2))',
                         repr(F))

    def test_lp_filter(self):
        F = operator.lp_filter(self.fine, self.coarse)
        repr(F)

    def test_injection_interpolation(self):
        P = operator.injection_interpolation(self.fine, self.coarse)
        Z = P.matching_zero()
        self.assertEqual('(interpolate\n  (grid 1 1)\n  (grid 2 2))', repr(P))
        self.assertEqual('0', repr(Z))

    def test_injection_restriction(self):
        R = operator.injection_restriction(self.fine, self.coarse)
        Z = R.matching_zero()
        self.assertEqual('(restrict\n  (grid 1 1)\n  (grid 2 2))', repr(R))
        self.assertEqual('0', repr(Z))

    def test_periodic_stencil(self):
        z = [ ((0,0), 0) ]
        flt = NdArray(dim=2, entries= [[ z ]])
        F = operator.from_periodic_stencil(flt, self.fine)
        self.assertEqual('(block\n  [[(stencil [((0, 0), 0.0)])]])', repr(F))

    def test_system(self):
        I = operator.identity(self.fine)
        Z = operator.zero(self.fine)
        L = system([[I, Z], [Z, I]])
        self.assertEqual('(system\n  [[id, 0], [0, id]])', repr(L))



