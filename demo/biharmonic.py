# Smoothing analysis of the biharmonic equation

from lfa_lab import *
import matplotlib.pyplot as mpp

grid = Grid(2)
coarse_grid = grid.coarse((2,2,))
Laplace = gallery.poisson_2d(grid)
I = operator.identity(grid)
Z = operator.zero(grid)

A = system([[Laplace, I],
            [Z      , Laplace]])
S_pointwise = jacobi(A, 0.8)
S_collective = collective_jacobi(A, 0.8)

symbol = S_pointwise.symbol()

print("Spectral radius: {}".format(symbol.spectral_radius()))
print("Spectral norm: {}".format(symbol.spectral_norm()))

plot.plot_2d(S_pointwise[0,0])
plot.plot_2d(S_pointwise[0,1])
plot.plot_2d(S_pointwise[1,0])
plot.plot_2d(S_pointwise[1,1])

plot.plot_2d(S_collective[0,0])
plot.plot_2d(S_collective[0,1])
plot.plot_2d(S_collective[1,0])
plot.plot_2d(S_collective[1,1])

# create diagonal filter
F = hp_filter(grid, coarse_grid)
FF = SystemNode([[F, Z], [Z, F]])

print((FF * S_pointwise).symbol().spectral_radius())
print((FF * S_collective).symbol().spectral_radius())

mpp.show()

