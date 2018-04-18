
import lfa_lab
from lfa_lab import *
from lfa_lab.plot import plot_2d
from lfa_lab.stencil import PeriodicStencil
import matplotlib.pyplot as plt

from lfa_lab.core import diff_jump_stride_2d, \
                       flux_conserving_int_2d
import lfa_lab.gallery as gallery

b = 4   # Block size
density1 = 1
density2 = 1e6

lfa_lab.dag.default_resolution = 128

assert(b % 2 == 0)

grid = Grid(2)
coarse_grid = grid.coarse((2,2))

# create the symbol
L_raw = diff_jump_stride_2d(grid.step_size(), b, density1, density2)
P_raw = flux_conserving_int_2d(L_raw)

def raw_pstencil_to_py(raw_pstencil, grid):
    """Convert a C++ periodic stencil into the Python version."""
    entries = NdArray(shape=raw_pstencil.shape())

    for p in NdRange(raw_pstencil.shape()):
        raw_stencil = raw_pstencil[p]
        stencil = SparseStencil()

        for q in ShiftedGrid(raw_stencil.startIndex(), raw_stencil.shape()):
            stencil.append(q, raw_stencil[q])

        entries[p] = stencil

    return PeriodicStencil(entries._entries)

st_L = raw_pstencil_to_py(L_raw, grid)
L = operator.from_periodic_stencil(st_L, grid)

# adaptive restriction/interpolation
st_P = raw_pstencil_to_py(P_raw, grid)
P_ad = operator.from_periodic_stencil(st_P, grid) \
    * FlatInterpolationNode(grid, coarse_grid)
R_ad = P_ad.adjoint()

# standard intergrid operators
P_st = gallery.ml_interpolation(grid, coarse_grid)
R_st = gallery.fw_restriction(grid, coarse_grid)

# RCA
#Lc_raw = diff_jump_stride_2d(coarse_grid.step_size(), b//2, 1, density2)
#st_Lc = raw_pstencil_to_py(Lc_raw, coarse_grid)
#Lc = st_Lc.symbol

# GCA
Lc_ad = R_ad * L * P_ad
Lc_st = R_st * L * P_st

S = jacobi(L, 0.8)
I = IdentityNode(grid)

print('Jacobi')
save_report(S, 'jc_smooth', 'JC Jacobi')

print('Standard CGC')
CGC_st = coarse_grid_correction(L, Lc_st, P_st, R_st)
save_report(CGC_st, 'jc_standard_cgc', 'JC Standard Coarse Grid Correction')

print('Adaptive CGC')
CGC_ad = coarse_grid_correction(L, Lc_ad, P_ad, R_ad)
save_report(CGC_ad, 'jc_adaptive_cgc', 'JC Adaptive Coarse Grid Correction')

# print(st_L)

print('Standard two-grid')
TG_st = S * CGC_st * S
save_report(TG_st, 'jc_standard_tg', 'JC Standard Two Grid Method')

print('Adaptive two-grid')
TG_ad = S * CGC_ad * S
save_report(TG_ad, 'jc_adaptive_tg', 'JC Adaptive Two Grid Method')

#plot_2d(L)
#plt.colorbar()
#plt.title('Operator')
print_report(L, 'Operator')

smpl = L.symbol()
print_report( (1.0 / smpl.norm()) * L * S, 'Hackbusch Smoothing')

#plot_2d(S)
#plt.colorbar()
#plt.title('Smoother')
print_report(S, title='JC Jacobi')

print_report(S * L.inverse(), 'JC Jacobi (res.)')

print_report(CGC_ad, 'Adaptive CGC')
print_report(CGC_st, 'Standard CGC')

#ews = S.symbol((32,32)).eigenvalues()
#plt.figure()
#plt.scatter(np.real(ews), np.imag(ews))

print_report(TG_ad, 'Adaptive TG')
print_report(TG_st, 'Standard TG')

plt.show()
