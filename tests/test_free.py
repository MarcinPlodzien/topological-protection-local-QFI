#!/usr/bin/env python3
"""Checks of codes/free_fermion.py:  python3 tests/test_free.py
 F1  window kernel vs explicit time integration
 F2  time-averaged QFI and QFI of the time-averaged state vs explicit time integration
 F3  clean gamma=1: stationary QFI == (1-h^2)^2, R_x stationary QFI == 0
 F4  encoding-site scan: stationary QFI(j=1,k) == (1-h^2)^2 h^(2(k-1))
 F7  single-spin QFI at read-out site j > 1 (Pfaffian, Jordan-Wigner string) vs exact diagonalization of the spin chain
 F5  Lyapunov criterion: clean chain exponent == ln(h/J) for gamma=1; sign change at h=J for gamma=0.5
 F6  disordered gamma=1 chain: closed form [sum_j prod_{l<j}(h_l/J)^2]^-2 == numerical stationary QFI (L=300, T=1e7);
     clean limit (1-h^2)^2 = 0.5625
Run from the tests folder or from the repository root; imports ../codes.  Runtime about 4 s (2 threads).
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes'))
from free_fermion import *
ok = True
def check(name, err, tol):
    global ok; good = err < tol; ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}: {err:.2e} (tol {tol:.0e})")
rng = np.random.default_rng(0)
L = 30; T = 40.0
h = 0.5 + rng.uniform(-1, 1, L)
w, v = spec(Amat(L, 1, 0.7, h))
ts = np.linspace(0, T, 8001)
tr = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
check("F1 kernel", abs(tr(np.exp(-1j * 0.37 * ts), ts) / T - window_kernel(0.37, T)), 1e-6)
Gs = np.array([Gt(w, v, t) for t in ts[::8]]); tt = ts[::8]
r = late_time_general(w, v, L, T)
rc = late_time(*chiral_modes(L, 1, 0.7, h), T)
for key in r: check(f'F2 chiral fast path == general ({key})', abs(r[key] - rc[key]), 1e-10)
rc3 = late_time(*chiral_modes(L, 1, 0.7, h), T, k=3); r3 = late_time_general(w, v, L, T, k=3)
for key in r3: check(f'F2 chiral fast path == general, k=4 ({key})', abs(r3[key] - rc3[key]), 1e-10)
for nm, src in (('y', 0), ('x', L)):
    avg = tr(Gs[:, 0, src] ** 2 + Gs[:, L, src] ** 2, tt) / T
    st = (tr(Gs[:, 0, src], tt) / T) ** 2 + (tr(Gs[:, L, src], tt) / T) ** 2
    check(f"F2 time-averaged QFI ({nm})", abs(avg - r['avg_' + nm]), 1e-5)
    check(f"F2 QFI of time-averaged state ({nm})", abs(st - r['stat_' + nm]), 1e-5)
for hh in (0.2, 0.5, 0.8):
    L = 200; M = chiral_modes(L, 1, 1, hh); r = late_time(*M, 2000.0)
    check(f"F3 stationary R_y, h={hh}", abs(r['stat_y'] - (1 - hh * hh) ** 2), 1e-4)
    check(f"F3 stationary R_x, h={hh}", r['stat_x'], 1e-6)
    for k in (1, 2, 4):
        check(f"F4 encoding site k={k+1}, h={hh}", abs(late_time(*M, 2000.0, k=k)['stat_y'] - (1 - hh * hh) ** 2 * hh ** (2 * k)), 1e-4)
check("F5 gamma=1 clean", abs(lyapunov_zero_mode(np.full(100, 0.5), 1, 1.0) - np.log(0.5)), 1e-12)
for hh, gm in ((0.9, 0.5), (1.1, 0.5), (0.3, 0.5), (1.5, 0.3)):      # exact roots of (1+g)x^2 - 2hx + (1-g) = 0
    ex = np.log(np.abs(np.roots([1 + gm, -2 * hh, 1 - gm])).max())
    check(f"F5 gamma={gm}, h={hh}: exponent {ex:+.3f} ({'topological' if ex < 0 else 'trivial'})", abs(lyapunov_zero_mode(np.full(20000, hh), 1, gm) - ex), 2e-3)
hd = 0.5 + np.random.default_rng(5).uniform(-1.2, 1.2, 300)
check("F6 disordered Ising chain: exact zero-mode formula vs numerics", abs(stationary_exact_ising(hd) - late_time(*chiral_modes(300, 1, 1.0, hd), 1e7)['stat_y']), 1e-6)
check("F6 clean limit", abs(stationary_exact_ising(np.full(4000, 0.5)) - 0.5625), 1e-12)
from ed_interacting import dense_traces
for (gm, hh, kk, tt_) in ((0.7, 0.5, 0, 3.0), (1.0, 0.3, 2, 5.5), (0.3, 1.4, 1, 2.0)):
    L7 = 8; w7, v7 = spec(Amat(L7, 1, gm, hh)); Fy, Fx = spin_qfi_profile(w7, v7, L7, tt_, k=kk)
    ed = np.array([dense_traces(L7, 1.0, gm, hh, 0.0, [tt_], k=kk, j=j)[0] for j in range(L7)])
    check(f"F7 spin QFI at all sites vs ED, gamma={gm}, h={hh}, k={kk + 1}", max(abs(ed[:, 0] - Fy).max(), abs(ed[:, 1] - Fx).max()), 1e-10)
print("\nALL TESTS PASSED" if ok else "\nSOME TESTS FAILED"); sys.exit(0 if ok else 1)
