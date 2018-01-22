from lfa_lab import *
import lfa_lab.gallery

grid = Grid(2)

stencil = lfa_lab.gallery.poisson_2d(grid)
J = jacobi(stencil, 0.8)

print_report(J)
save_report(J, 'test_report', 'Jacobi Smoother on Poisson', True)
