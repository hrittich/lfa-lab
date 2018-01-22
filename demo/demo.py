# -*- coding: utf-8 -*-
from lfa_lab import *
import matplotlib.pyplot as plt

plot.default_options['style_2d'] = 'mesh'

# Create a 2d grid with default spacing
grid = Grid(2)
coarse = grid.coarse((2,2))

A = gallery.poisson_2d(grid)

plot.plot_2d(A)

E = jacobi(A, 0.8)
plot.plot_2d(E)

S = HpFilterNode(grid, coarse);
plot.plot_2d( S * E )

plt.show()
