#!/usr/bin/env python3
"""
run_spin_maps.py -- exact single-SPIN QFI maps F_j(t) (Jordan-Wigner string included, Pfaffian formula)
for FigR1(b,c), plus the fermionic-mode expression |U+V|^2 for comparison (they coincide only at j=1).
Writes data/free/spinmap_{tag}.txt (rows: time, columns: site) and data/free/spinmap_times.txt.
    python3 codes/run_spin_maps.py [nproc]

Parameters: L = 100, rotation at site 1 (k = 0 in the code), read-out sites j = 1..60, 121 times tJ = linspace(0, 60, 121);
tag = topological (gamma = 0.3, h = 0.5J) and trivial (gamma = 0.8, h = 1.8J); J = 1, delta = 0, all spins down.
Outputs (data/free/):
  spinmap_times.txt               the 121 times
  spinmap_{tag}_Ry.txt            121 x 60, exact single-spin QFI for R_y (FigR1b,c; FigR2)
  spinmap_{tag}_Rx.txt            121 x 60, the same for R_x (not plotted)
  spinmap_{tag}_fermionic.txt     121 x 60, |U_j^1 + V_j^1|^2 = G[a_j,a_1]^2 + G[b_j,a_1]^2 (FigR2)
For read-out site j (1-based) the Pfaffian of a 2j x 2j matrix (up to 120 x 120) is evaluated for every time and axis
(free_fermion.spin_qfi_profile).
Multiprocessing pool of nproc processes (default 4), one thread each.  Deterministic.
Runtime 95 s with 2 processes, peak memory of the parent 32 MB.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from free_fermion import *
OUT = os.path.join(HERE, '..', 'data', 'free'); os.makedirs(OUT, exist_ok=True)
L = 100; JMAX = 60; TS = np.linspace(0, 60, 121)
CASES = {'topological': (0.3, 0.5), 'trivial': (0.8, 1.8)}

def snap(args):
    tag, t = args; g, h = CASES[tag]
    w, v = spec(Amat(L, 1.0, g, h))
    Fy, Fx = spin_qfi_profile(w, v, L, t, k=0, jmax=JMAX)
    G = Gt(w, v, t)
    return Fy, Fx, G[:JMAX, 0] ** 2 + G[L:L + JMAX, 0] ** 2

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    np.savetxt(os.path.join(OUT, "spinmap_times.txt"), TS, header="time grid tJ")
    with Pool(nproc) as p:
        for tag, (g, h) in CASES.items():
            r = p.map(snap, [(tag, t) for t in TS], chunksize=1)
            hd = f"L={L} gamma={g} h={h}, rotation at site 1; rows: times (spinmap_times.txt), columns: read-out site j=1..{JMAX}\n"
            np.savetxt(os.path.join(OUT, f"spinmap_{tag}_Ry.txt"), [x[0] for x in r], fmt="%.6e", header=hd + "exact single-spin QFI, R_y")
            np.savetxt(os.path.join(OUT, f"spinmap_{tag}_Rx.txt"), [x[1] for x in r], fmt="%.6e", header=hd + "exact single-spin QFI, R_x")
            np.savetxt(os.path.join(OUT, f"spinmap_{tag}_fermionic.txt"), [x[2] for x in r], fmt="%.6e", header=hd + "|U+V|^2 (fermionic-mode expression; equals the spin QFI only at j=1)")
            print("wrote", tag, flush=True)
