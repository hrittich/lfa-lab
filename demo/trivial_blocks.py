"""Trivial block demo.

Should produce the same plot multiple times.
"""
import lfa_lab
import lfa_lab.plot
import lfa_lab.gallery
from lfa_lab import NdArray

import matplotlib.pyplot as plt

grid = lfa_lab.Grid(2)
L = lfa_lab.gallery.poisson_2d(grid)
lfa_lab.plot.plot_2d(L)
plt.title('FD Poisson')


L2 = lfa_lab.BlockNode( [[L, L], \
                     [L, L]] )
lfa_lab.plot.plot_2d(L2)
plt.title('Poisson Block Test')

stencil = L.stencil
L3 = lfa_lab.operator.from_periodic_stencil(
    NdArray(dim=2, entries=
        [[stencil, stencil],
        [stencil, stencil]]),
    grid)
lfa_lab.plot.plot_2d(L3)
plt.title('Poisson Block Stencil Test')

plt.show()

