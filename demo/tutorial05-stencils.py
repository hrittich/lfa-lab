from lfa_lab import *

fine = Grid(2, [1.0/32, 1.0/32])
coarse = fine.coarse((2,2))

# A simple function that computes the poisson stencil for a given grid.
def my_poisson_2d(grid):
    h1, h2 = grid.step_size()
    entries = [
        (( 0, -1), -1.0 / (h2*h2)),
        ((-1,  0), -1.0 / (h1*h1)),
        (( 0,  0),  2.0 / (h1*h1) + 2.0 / (h2*h2)),
        (( 1,  0), -1.0 / (h1*h1)),
        (( 0,  1), -1.0 / (h2*h2))
      ]
    return operator.from_stencil(entries, grid)

L = my_poisson_2d(fine)
Lc = my_poisson_2d(coarse)
S = jacobi(L, 0.8)

# Create the restriction operator
entries = [
    ((-1, -1), 1.0/16),
    (( 0, -1), 1.0/8),
    (( 1, -1), 1.0/16),
    ((-1,  0), 1.0/8),
    (( 0,  0), 1.0/4),
    (( 1,  0), 1.0/8),
    ((-1,  1), 1.0/16),
    (( 0,  1), 1.0/8),
    (( 1,  1), 1.0/16),
]
R = operator.injection_restriction(fine, coarse) * \
        operator.from_stencil(entries, fine)

# Create the interpolation operator.
entries = [
    ((-1, -1), 1.0/4),
    (( 0, -1), 1.0/2),
    (( 1, -1), 1.0/4),
    ((-1,  0), 1.0/2),
    (( 0,  0), 1.0),
    (( 1,  0), 1.0/2),
    ((-1,  1), 1.0/4),
    (( 0,  1), 1.0/2),
    (( 1,  1), 1.0/4),
  ]
P = operator.from_stencil(entries, fine) * \
        operator.injection_interpolation(fine, coarse)

cgc = coarse_grid_correction(
        operator = L,
        coarse_operator = Lc,
        interpolation = P,
        restriction = R)
E = S * cgc * S

print(E.symbol().spectral_radius())
