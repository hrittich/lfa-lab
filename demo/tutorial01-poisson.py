from lfa_lab import *

# Create a 2D grid with step-size (1/32, 1/32).
fine = Grid(2, [1.0/32, 1.0/32])

# Create a poisson operator.
L = gallery.poisson_2d(fine)

# Create the jacobi smoother for L.
S = jacobi(L, 0.8)

# Compute the symbol of S.
symbol = S.symbol()

# Print spectral radius of S.
print(symbol.spectral_radius())
