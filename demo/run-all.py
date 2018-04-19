import matplotlib

# non interactive graphics
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# turn interactive mode on, such that plt.show will not block
plt.ion()

scripts = [
    'demo.py',
    'tutorial01-poisson.py',
    'tutorial02-plot.py',
    'tutorial03-smoothing.py',
    'tutorial04-twogrid.py',
    'tutorial05-stencils.py',
    'tutorial06-multigrid.py',
    'aggressive_coarse.py',
#    'biharmonic.py',
    'block_jac_2x2.py',
    'block_jac.py',
    'block_rb_jac_2x2.py',
    'jump_coeff_2d.py',
    'poisson.py',
    'rb_jacobi.py',
#    'red_black_block.py',
    'report.py',
    'trivial_blocks.py',
#    'unit-tests.py',
#    'var_coeff.py',
]

for s in scripts:
    print('=' * 80)
    print('Running {}'.format(s))
    with open(s, 'r') as fp:
        exec(fp.read())

    plt.close('all')


