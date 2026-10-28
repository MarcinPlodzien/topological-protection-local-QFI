"""Edge cases of the one-axis diagnostic (reply).
(1) Trivial chain (h_j = 1.8J) with one vanishing field h_m = 0: exact zero mode and one-axis response at gamma = 1 only.
(2) Class-D term i*kappa*a_1*a_2 that breaks the chiral symmetry: the stationary signal survives, the one-axis property does not.
Run: python3 codes/check_edge_cases.py > data/free/edge_cases.txt

Model: free XY chain (delta = 0), J = 1, L = 200, k = j = 1, all spins down; stationary QFI from free_fermion.late_time
(window [0, T]).  Part (2) adds i*kappa*a_1*a_2 to the Majorana matrix of bdg.Amat (A[0,1] += 2 kappa), which breaks
the chiral structure, so the 2L x 2L reference late_time_general is used there.
Output (stdout, two blocks with their own comment line):
  (1) gamma m h_m lowest|eps| then stat_y stat_x at TJ = 1e2, 1e4, 1e6   (gamma = 1, 0.8, 0.5; h_2 = 0, h_3 = 0, h_2 = 0.05)
  (2) kappa lowest|eps| then stat_y stat_x at TJ = 1e3, 1e5              (kappa = 0, 0.1, 0.5)
Quoted in the reply.  Deterministic.  Runtime 0.8 s (2 threads).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from free_fermion import chiral_modes, late_time, late_time_general
from bdg import Amat, spec

L = 200
print(f"# (1) L={L}, h_j=1.8J except h_m; columns: gamma, m, h_m, lowest |eps|, then stat_y stat_x for TJ=1e2, 1e4, 1e6")
for gam in (1.0, 0.8, 0.5):
    for m, hm in ((2, 0.0), (3, 0.0), (2, 0.05)):
        if gam != 1.0 and hm != 0.0:
            continue
        h = np.full(L, 1.8); h[m - 1] = hm
        U, s, V = chiral_modes(L, 1.0, gam, h)
        row = [late_time(U, s, V, T) for T in (1e2, 1e4, 1e6)]
        print(f"{gam:.1f} {m} {hm:.2f} {s.min():.3e} " + " ".join(f"{r['stat_y']:.4e} {r['stat_x']:.4e}" for r in row))
print(f"# (2) L={L}, gamma=1, h=0.5J, term i*kappa*a_1*a_2; columns: kappa, lowest |eps|, then stat_y stat_x for TJ=1e3, 1e5")
for kap in (0.0, 0.1, 0.5):
    A = Amat(L, 1.0, 1.0, 0.5); A[0, 1] += 2 * kap; A[1, 0] -= 2 * kap
    w, v = spec(A)
    row = [late_time_general(w, v, L, T) for T in (1e3, 1e5)]
    print(f"{kap:.1f} {np.abs(w).min():.3e} " + " ".join(f"{r['stat_y']:.4e} {r['stat_x']:.4e}" for r in row))
