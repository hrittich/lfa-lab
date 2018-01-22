from lfa_lab import *

def multigrid(level, fine_grid):
    """Return a tuple consisting of the error propagation operator and the
    linear system operator of a multigrid method."""
    L = gallery.poisson_2d(fine_grid)

    if level == 1:
        # Solve exactly on the coarsest grid, hence return the zero operator.
        E = operator.zero(fine_grid)
        return (E, L)
    else:
        coarse_grid = fine_grid.coarse((2,2))
        Ec, Lc = multigrid(level-1, coarse_grid)

        S = jacobi(L, 0.8)

        R = gallery.fw_restriction(fine_grid, coarse_grid)
        P = gallery.ml_interpolation(fine_grid, coarse_grid)

        cgc = coarse_grid_correction(
                operator = L,
                coarse_operator = Lc,
                interpolation = P,
                restriction = R,
                coarse_error = Ec)

        # Apply one pre- and one post-smoothing step.
        E = S * cgc * S

        return (E, L)

# Compute error propagation operator of a three-grid method
fine_grid = Grid(2, [1.0/32, 1.0/32])
E, L = multigrid(3, fine_grid)

print(E.symbol().spectral_radius())
