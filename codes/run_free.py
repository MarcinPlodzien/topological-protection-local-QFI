#!/usr/bin/env python3
"""
run_free.py -- all clean and engineered free-fermion data of the paper (delta = 0).

Model: H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - sum h_j Z (Pauli), OBC, J = 1, initial state all down,
rotation R_y or R_x at site k (k = 1 unless stated), read-out of the edge spin j = 1, theta = 0.
"stat" = QFI of the state averaged over the window [0, T] (stationary QFI), "avg" = QFI averaged over [0, T];
both from the closed-form window kernel of free_fermion.late_time (no time grid).

Outputs (data/free/, plain text with a header line):
  A. window_scan_L{8,12,16,20,200}.txt   gamma=1, h=0.5J; T = logspace(0, 9, 73)
                                         columns: T stat_y stat_x avg_y avg_x          (FigR4b)
     realtime_L400_g1.00_h0.50.txt       L=400, no averaging, tJ = linspace(0, 60, 1201): tJ F_y F_x  (FigR4a)
     realtime_L400_g1.00_h1.80.txt       the same for h=1.8J (trivial), FigR4a inset
  B. encoding_scan.txt                   L=200, TJ=1e3, (gamma,h) = (1,0.3), (1,0.6), (0.5,0.6), k = 1..12
                                         columns: gamma h k stat_y stat_x avg_y avg_x exact(gamma=1)  (FigR4c)
  C. plateau_vs_h.txt                    L=400, TJ=2e3, gamma = 1, 0.5, 0.25, h = linspace(0, 1.4, 71)
                                         columns: gamma h stat_y stat_x avg_y avg_x exact_semi_infinite  (FigR4d)
  D. cases_{tag}.txt, spectrum_{tag}.txt L=200, gamma=1, T = logspace(0.5, 10, 77)
                                         tag = topological (h=0.5J), trivial_clean (1.8J), trivial_impurity (h_1=0.3J),
                                         trivial_finetuned (h_1=0), island_n{4,8,16} (first n sites 0.5J, rest 1.8J),
                                         smooth_ramp (linspace(0.3, 1.8, 30) on the first 30 sites, then 1.8J)
                                         cases: T stat_y stat_x avg_y avg_x; header: three lowest single-particle energies
                                         spectrum: s_n, U_1n^2, V_1n U_1n, V_1n^2 (sorted by s_n)  (Figs. 3b, S2b, TableR1)
     island_splitting.txt                n = 2..30: n eps0 stat_y(T=30) stat_x(T=30)  (FigR8a)

Run (from the repository root):  python3 codes/run_free.py
Single process; sets OMP/OpenBLAS threads to 2 unless already set.  Deterministic (no random numbers).
Runtime 292 s, peak memory 98 MB (2 threads).
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from free_fermion import *
OUT = os.path.join(HERE, '..', 'data', 'free'); os.makedirs(OUT, exist_ok=True)
J = 1.0
def save(name, arr, header):
    np.savetxt(os.path.join(OUT, name), arr, fmt="%.8e", header=header); print("wrote", name, flush=True)

# --- A. stationary state and order of limits: window [0,T], gamma=1, h=0.5 -----------------------
Ts = np.logspace(0, 9, 73)
for L in (8, 12, 16, 20, 200):
    M = chiral_modes(L, J, 1.0, 0.5); e0 = M[1].min()
    rows = [[T] + [late_time(*M, T)[k] for k in ('stat_y', 'stat_x', 'avg_y', 'avg_x')] for T in Ts]
    save(f"window_scan_L{L}.txt", rows, f"gamma=1 h=0.5 L={L}, edge splitting eps0={e0:.6e}; window [0,T]\nT  stat_y stat_x avg_y avg_x")
# real-time trace, no averaging, thermodynamic limit (L=400, t < L/v)
L = 400; w, v = spec(Amat(L, J, 1.0, 0.5)); ts = np.linspace(0, 60, 1201)
save("realtime_L400_g1.00_h0.50.txt", [[t, *qfi_t(w, v, L, t)] for t in ts], "gamma=1 h=0.5 L=400, no averaging\ntJ F_y F_x")
w, v = spec(Amat(L, J, 1.0, 1.8))                                   # the same in the trivial phase
save("realtime_L400_g1.00_h1.80.txt", [[t, *qfi_t(w, v, L, t)] for t in ts], "gamma=1 h=1.8 L=400, no averaging\ntJ F_y F_x")

# --- B. encoding-site scan: read-out j=1, rotation at site k ------------------------------------
L = 200; T = 1e3; rows = []
for (g, h) in ((1.0, 0.3), (1.0, 0.6), (0.5, 0.6)):
    M = chiral_modes(L, J, g, h)
    for k in range(12):
        r = late_time(*M, T, k=k)
        rows.append([g, h, k + 1, r['stat_y'], r['stat_x'], r['avg_y'], r['avg_x'], (1 - h * h) ** 2 * h ** (2 * k) if g == 1 else np.nan])
save("encoding_scan.txt", rows, "L=200, read-out site 1, window [0,T], T=1e3\ngamma h k stat_y stat_x avg_y avg_x exact(gamma=1)=(1-h^2)^2 h^(2(k-1))")

# --- C. plateau versus field for several anisotropies ---------------------------------------------
hs = np.linspace(0.0, 1.4, 71); rows = []
for g in (1.0, 0.5, 0.25):
    for h in hs:
        L = 400; r = late_time(*chiral_modes(L, J, g, h), 2e3)
        rows.append([g, h, r['stat_y'], r['stat_x'], r['avg_y'], r['avg_x'], plateau_exact(h, J, g)])
save("plateau_vs_h.txt", rows, "L=400, k=j=1, window [0,T], T=2e3\ngamma h stat_y stat_x avg_y avg_x exact_semi_infinite")

# --- D. topological edge mode versus trivially localized states -------------------------------
L = 200; Ts = np.logspace(0.5, 10, 77)
def case(tag, h, g=1.0):
    M = chiral_modes(L, J, g, h); U, s_, V = M; aw = np.sort(s_)[:3]
    rows = [[T] + [late_time(*M, T)[k] for k in ('stat_y', 'stat_x', 'avg_y', 'avg_x')] for T in Ts]
    save(f"cases_{tag}.txt", rows, f"L=200 gamma={g}; three lowest |omega|: {aw[0]:.4e} {aw[1]:.4e} {aw[2]:.4e}\nT stat_y stat_x avg_y avg_x")
    # edge response at site 1: F_y(t) = (sum_n U_1n^2 cos s_n t)^2 + (sum_n V_1n U_1n sin s_n t)^2
    o = np.argsort(s_)
    save(f"spectrum_{tag}.txt", np.column_stack([s_[o], U[0, o] ** 2, V[0, o] * U[0, o], V[0, o] ** 2]),
         "mode frequency s_n | U_1n^2 (R_y, cos) | V_1n U_1n (sin, both axes) | V_1n^2 (R_x, cos)")
case("topological", np.full(L, 0.5))
h = np.full(L, 1.8); case("trivial_clean", h)
h = np.full(L, 1.8); h[0] = 0.3; case("trivial_impurity", h)
h = np.full(L, 1.8); h[0] = 0.0; case("trivial_finetuned", h)
for n in (4, 8, 16):
    h = np.full(L, 1.8); h[:n] = 0.5; case(f"island_n{n}", h)
h = np.full(L, 1.8); h[:30] = np.linspace(0.3, 1.8, 30); case("smooth_ramp", h)
# island: lifetime 1/eps0 versus island length
rows = []
for n in range(2, 31):
    h = np.full(L, 1.8); h[:n] = 0.5; M = chiral_modes(L, J, 1.0, h); r = late_time(*M, 30.0); rows.append([n, M[1].min(), r['stat_y'], r['stat_x']])
save("island_splitting.txt", rows, "first n sites h=0.5 (topological), rest h=1.8 (trivial), gamma=1, L=200\nn eps0 stat_y(T=30) stat_x(T=30)")
