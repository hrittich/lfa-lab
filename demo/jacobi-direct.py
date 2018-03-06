from lfa_lab import *
import matplotlib.pyplot as mpp

grid = Grid(2, [1.0/32, 1.0/32])
A = gallery.poisson_2d(grid)
I = operator.identity(grid)
omega = 0.8
E = I - omega * A.diag().inverse() * A

plot.plot_2d(E.symbol())
mpp.show()
