"""Trivial block demo.

Should produce the same plot multiple times.
"""
import lfa_lab
import lfa_lab.plot
import lfa_lab.gallery
from lfa_lab import NdArray

import matplotlib.pyplot as plt

# We compute the well-known Poisson symbol
grid = lfa_lab.Grid(2)
L = lfa_lab.gallery.poisson_2d(grid)
lfa_lab.plot.plot_2d(L)
plt.title('FD Poisson')

# Given a square pattern of operators, the BlockNode repeats this pattern
# over the whole domain and applies it to the grid. Since we are creating
# a pattern with the same operator everywhere, we should get the same result
# as before.
L2 = lfa_lab.BlockNode([[L, L], \
                        [L, L]])
lfa_lab.plot.plot_2d(L2)
plt.title('Poisson Block Test')

# Essentially the same thing as above, but constructing the operator directly
# from the stencils.
stencil = L.stencil
L3 = lfa_lab.operator.from_periodic_stencil(
    NdArray(dim=2,
            entries=[[stencil, stencil],
                     [stencil, stencil]]),
    grid)
lfa_lab.plot.plot_2d(L3)
plt.title('Poisson Block Stencil Test')

plt.show()

