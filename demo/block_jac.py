
from lfa_lab import *
from lfa_lab.plot import *
from lfa_lab.gallery import *

import matplotlib.pyplot as plt

fine = Grid(2)
coarse = fine.coarse((2,2))

L_st = poisson_2d(fine)

F = HpFilterNode(fine, coarse)
E = block_jacobi(L_st, (4,4), 0.7)

print_report(F*E)
print_report(E)

plt.show()
