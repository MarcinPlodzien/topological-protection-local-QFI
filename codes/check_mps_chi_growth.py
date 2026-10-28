#!/usr/bin/env python3
"""
check_mps_chi_growth.py -- growth of the bond dimension in the matrix-product-state runs (reply,
Sec. VII, "Long chains"; reply FigR12).

Reads data/mps/mps_L48_g1.00_h0.20_d0.20_chi128_dt0.1.txt and mps_L100_..._d0.20_... (columns tJ F_y F_x chi_psi chi_phi
F_global norm_phi), fits ln chi_phi = a + r t on 2 <= tJ <= 6 (before the cap chi = 128), and prints the rate r,
the doubling time ln2/r and the times at which the fitted growth would reach chi = 10^3 and 10^4 (an extrapolation).
Also prints max |F_y(L=100) - F_y(L=48)| for tJ <= 8 and whether chi is equal at all common times.
Run: python3 codes/check_mps_chi_growth.py > data/mps/chi_growth.txt   (seconds)
"""
import os
import numpy as np

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'mps')
a = np.loadtxt(os.path.join(D, 'mps_L48_g1.00_h0.20_d0.20_chi128_dt0.1.txt'))
b = np.loadtxt(os.path.join(D, 'mps_L100_g1.00_h0.20_d0.20_chi128_dt0.1.txt'))
for name, m in (('L=48', a), ('L=100', b)):
    k = (m[:, 0] >= 2) & (m[:, 0] <= 6); r, c = np.polyfit(m[k, 0], np.log(m[k, 4]), 1)
    print(f"{name}: rate {r:.3f}/J, doubling time {np.log(2) / r:.2f}/J, chi=1e3 at tJ {(np.log(1e3) - c) / r:.1f}, "
          f"chi=1e4 at tJ {(np.log(1e4) - c) / r:.1f}")
n = min(len(a), len(b)); same = np.array_equal(a[:n, 4], b[:n, 4]) and np.allclose(a[:n, 0], b[:n, 0])
k = a[:n, 0] <= 8
print(f"chi equal at all {n} common times: {same}; max |F_y(100)-F_y(48)| for tJ<=8: {np.abs(a[:n, 1] - b[:n, 1])[k].max():.2e}")
