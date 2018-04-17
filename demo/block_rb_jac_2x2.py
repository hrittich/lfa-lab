from lfa_lab import *
import matplotlib.pyplot as plt

grid = Grid(2, [1.0/32, 1.0/32])

# Define the discrete Laplace operator
a = [ ((-1,0),-1), ((0, -1), -1),
      ((0,0), 4), ((0,1), -1), ((1,0), -1) ]
A = operator.from_stencil(a, grid)

# Define the block diagonal part of A
d = NdArray(shape=(2,2))
d[0,0] = [ ((0,0), 4), ((0,1), -1), ((1,0), -1) ]
d[0,1] = [ ((0,-1), -1), ((0,0), 4), ((1,0), -1) ]
d[1,0] = [ ((-1,0), -1), ((0,0), 4), ((0,1), -1) ]
d[1,1] = [ ((-1,0), -1), ((0,-1), -1), ((0,0), 4) ]
D = operator.from_periodic_stencil(d, grid)

# Define the filtering operators
z = [ ((0,0), 0) ]
o = [ ((0,0), 1) ]
red = NdArray(dim=2, entries=
    [[ o, o, z, z ],
     [ o, o, z, z ],
     [ z, z, o, o ],
     [ z, z, o, o ]])
FR = from_periodic_stencil(red, grid)
black = NdArray(dim=2, entries=
    [[ z, z, o, o ],
     [ z, z, o, o ],
     [ o, o, z, z ],
     [ o, o, z, z ]])
FB = from_periodic_stencil(black, grid)

# Define the block Jacobi error propagation operator
I = operator.identity(grid)
J = (I - D.inverse() * A)

# Define the RB-Jacobi error propagator
E = (FR + FB * J) * (FB + FR * J)

print((E.symbol().spectral_radius()))

plot.plot_2d(E, norm_type='output')
plt.show()




