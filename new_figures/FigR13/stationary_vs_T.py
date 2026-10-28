"""Stationary QFI of the edge spin versus averaging time T, full diagonalization, L = 12, gamma = 1.
Writes stationary_vs_T_<case>.txt with columns T, I_st(R_y), I_st(R_x)."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stationary_qfi_ed import stationary
L = 12; Ts = np.logspace(0, 4, 21)
cases = {"topological_h0.20_d0.10": (np.full(L, 0.2), 0.1),
         "trivial_h1_0.00_d0.00": (np.r_[0.0, np.full(L - 1, 1.8)], 0.0),
         "trivial_h1_0.00_d0.10": (np.r_[0.0, np.full(L - 1, 1.8)], 0.1)}
for name, (h, d) in cases.items():
    r = stationary(L, 1.0, h, d, Ts)
    np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"stationary_vs_T_{name}.txt"), [[T, r[T][0], r[T][1]] for T in Ts], fmt="%.6e",
               header=f"L={L} gamma=1 fields={np.round(h,3).tolist()} delta={d}; stationary QFI of the edge spin\nT I_st_y I_st_x")
    print(name, "done", flush=True)
