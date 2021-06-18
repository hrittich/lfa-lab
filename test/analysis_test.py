from unittest import main, TestCase

import numpy as np
from lfa_lab import *

class AnalysisTest(TestCase):

    def test_poisson_h_ellipticity(self):
        fine = Grid(2, [1.0/32, 1.0/32])
        L = gallery.poisson_2d(fine)

        self.assertAlmostEqual(h_ellipticity(L), 0.25)

    def test_poisson_smoothing_factor(self):
        fine = Grid(2, [1.0/32, 1.0/32])
        L = gallery.poisson_2d(fine)
        J = smoother.jacobi(L, 4.0/5.0)

        self.assertLess(abs(smoothing_factor(J) - 3.0/5), 1e-2)


if __name__ == '__main__':
    main()
