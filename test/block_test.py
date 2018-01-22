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
#import lfa_lab.plot

import unittest
from unittest import TestCase

import numpy as np

class Block(TestCase):

    def test_injection(self):

        grid = Grid(2)

        ops = NdArray(shape=(4,4))

        for i in NdRange(ops.shape):
            i = np.array(i)
            if (i == (1,0)).all():
                ops[i] = IdentityNode(grid)
            else:
                ops[i] = ZeroNode(grid)

        op = BlockNode(ops._entries)
        sym = op.symbol( (4, 4) )
        #print(sym)
        #lfa_lab.plot.plot_symbol_2d(sym)
        #lfa_lab.plot.plot_fo_2d(op)

        # This test is incomplete.
        # TODO
        # Check result...


if __name__ == '__main__':
    unittest.main()

