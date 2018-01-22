from lfa_lab import *
from lfa_lab.plot import plot_2d
import lfa_lab.gallery as gallery
import matplotlib.pyplot as plt

def ac_two_grid(coarsening, block_size):
    """Aggressive coarsening two grid."""

    grid = Grid(2)
    coarse_grid = grid.coarse(coarsening)

    A = gallery.poisson_2d(grid)
    Ac = gallery.poisson_2d(coarse_grid)

    R = gallery.fw_restriction(grid, coarse_grid)
    P = gallery.ml_interpolation(grid, coarse_grid)

    cgc = coarse_grid_correction(
            operator = A,
            coarse_operator = Ac,
            interpolation = P,
            restriction = R)

    J = block_jacobi(A, block_size, .8)
    # Two grid method
    TG = J * cgc * J
    return TG

plot_2d(ac_two_grid((2,2), (1,1)))
plt.show()

