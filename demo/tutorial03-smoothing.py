from lfa_lab import *

fine = Grid(2, [1.0/32, 1.0/32])

# Create a coarse grid by selecting every second point in the x- and every
# second point in the y-direction.
coarse = fine.coarse((2,2))

L = gallery.poisson_2d(fine)
S = jacobi(L, 0.8)

# Create the filtering operator.
F = operator.hp_filter(fine, coarse)

# Apply the filter to the smoother. After applying smoother all low modes are
# removed from the spectrum.
E = F * S

print(E.symbol().spectral_radius())
