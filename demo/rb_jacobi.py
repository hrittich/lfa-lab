from lfa_lab import *
import matplotlib.pyplot as plt

g = Grid(2)
s = gallery.poisson_2d(g)
S = rb_jacobi(s)

print_report(S)
plt.show()

