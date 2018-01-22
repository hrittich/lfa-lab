from lfa_lab import *
# To plot a symbol, we need to import matplotlib
import matplotlib.pyplot as mpp

fine = Grid(2, [1.0/32, 1.0/32])
L = gallery.poisson_2d(fine)
S = jacobi(L, 0.8)
symbol = S.symbol()

# Plot the symbol
plot.plot_2d(symbol)

# Tell matplotlib to show the plot.
mpp.show()
