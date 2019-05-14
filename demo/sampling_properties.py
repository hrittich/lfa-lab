import lfa_lab
import lfa_lab.gallery
from math import pi
import numpy as np

grid = lfa_lab.Grid(2)

L = lfa_lab.gallery.poisson_2d(grid)
S = lfa_lab.jacobi(L, .8)

# Sample with a zero base frequency
symbol = S.symbol(base_frequency = (0,0))
print(symbol.spectral_radius())

# Sample with a base frequency that is as far from zero as possible.
# (This is the default when no base frequency is given.)
desired_resolution = (128, 128)
resolution = S.properties.adjustResolution(desired_resolution)
base_frequency = pi / (np.array(grid.step_size())
                       * np.array(resolution))
symbol = S.symbol(desired_resolution=resolution, base_frequency=base_frequency)
print(symbol.spectral_radius())


