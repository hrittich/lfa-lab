import lfa_lab
import lfa_lab.gallery
import lfa_lab.plot
import matplotlib.pyplot as plt

grid = lfa_lab.Grid(2)
coarse_grid = grid.coarse((2,2))

L = lfa_lab.gallery.poisson_2d(grid)
Lc = lfa_lab.gallery.poisson_2d(coarse_grid)

S = lfa_lab.jacobi(L, .8)
R = lfa_lab.gallery.fw_restriction(grid, coarse_grid)
P = lfa_lab.gallery.ml_interpolation(grid, coarse_grid)

cgc = lfa_lab.coarse_grid_correction(
        operator = L,
        coarse_operator = Lc,
        interpolation = P,
        restriction = R)

E = S * cgc * S

#print('Input {0}'.format(E.inputModes()))
#print('Output {0}'.format(E.outputModes()))

lfa_lab.print_report(E)


# 1D
grid = lfa_lab.Grid(1)
coarse_grid = grid.coarse((2,))

L = lfa_lab.gallery.poisson_1d(grid)
lfa_lab.print_report(L)

plt.show()

