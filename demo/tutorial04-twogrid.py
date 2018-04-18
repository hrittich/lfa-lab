from lfa_lab import *

fine = Grid(2, [1.0/32, 1.0/32])
coarse = fine.coarse((2,2))

L = gallery.poisson_2d(fine)
Lc = gallery.poisson_2d(coarse)
S = jacobi(L, 0.8)

# Create restriction and interpolation operators.
R = gallery.fw_restriction(fine, coarse)
P = gallery.ml_interpolation(fine, coarse)

# Construct the coarse grid correction from the individual operators.
cgc = coarse_grid_correction(
        operator = L,
        coarse_operator = Lc,
        interpolation = P,
        restriction = R)

# Apply one pre- and one post-smoothing step.
E = S * cgc * S

print((E.symbol().spectral_radius()))
