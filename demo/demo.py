from lfa_lab import *
import matplotlib.pyplot as plt

plot.default_options['style_2d'] = 'mesh'

# Create a 2d grid with default spacing
grid = Grid(2)
coarse = grid.coarse((2,2))

# The gallery module contains stencils from applications. We create the
# stencil that stems from the finite-difference discretization of the
# Laplace operator.
A = gallery.poisson_2d(grid)

# We can inspect how the operator acts on wave functions of different
# frequency using the plot function.
plot.plot_2d(A)

# From the operator A, we can create the iteration operator of the Jacobi
# method that would solve a linear system involving A.
E = jacobi(A, 0.8)

# Note that when using LFA Lab, we actually define the formula for which we
# want to compute the symbol. This formula can be inspected by printing
# the operator of usint the repr function.
print(E)

# We inspect the iteration operator.
plot.plot_2d(E)

# A smoothing analysis assumes that we only have to consider the response of
# the operator to high frequencies. Hence, using a high-pass filter, we can
# conduct a smoothing analysis.
S = HpFilterNode(grid, coarse);
plot.plot_2d( S * E )

plt.show()

