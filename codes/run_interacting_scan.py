#!/usr/bin/env python3
"""
run_interacting_scan.py -- fine scan of the interaction strength with full diagonalization (L = 10, 12).

For each delta the edge QFI (R_y and R_x, k = j = 1) is averaged over the windows listed in WINDOWS.
Purpose: resolve the dependence of the plateau lifetime on delta, including resonances between the
field and the interaction, and the integrable line h = 0.
Writes data/interacting/scan_L{L}_g{gamma}_h{h}.txt:  delta | <F_y> per window | <F_x> per window
    python3 codes/run_interacting_scan.py

Model: as in run_interacting.py (Pauli, OBC, J = 1, all spins down, k = j = 1).  Full diagonalization
(numpy eigh of the dense 2^L x 2^L matrix); for each delta in -0.6 ... 0.6 (step 0.025, 49 values) and each window,
F_y and F_x are averaged over NPT = 120 times drawn uniformly in the window from a new np.random.default_rng(1)
for every window, so the sampling times are the same for every delta.
Scans: (L, gamma, h) = (10, 1, 0.2), (12, 1, 0.2), (12, 1, 0.4), (12, 1, 0), (10, 1, 0.4), (12, 0.5, 0.2).
Output columns: delta | <F_y> over the six windows | <F_x> over the six windows (13 columns).
FigR10(c) uses scan_L10_g1.00_h0.20.txt and scan_L12_g1.00_h0.20.txt, columns 1-4 (windows up to 1e3/J).
Runtime: about 10 h for all six scans (dominated by L = 12); 2 threads.
Memory: about 1 GB at L = 12 (dense 4096 x 4096 matrices).
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed_interacting as ed
OUT = os.path.join(HERE, '..', 'data', 'interacting'); os.makedirs(OUT, exist_ok=True)
WINDOWS = [(10, 30), (30, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 10000)]
NPT = 120
DELTAS = np.round(np.arange(-0.6, 0.6001, 0.025), 3)

def scan(L, g, h):
    f = os.path.join(OUT, f"scan_L{L}_g{g:.2f}_h{h:.2f}.txt")
    if os.path.exists(f): return
    rows = []
    for d in DELTAS:
        diag, off = ed.build_H(L, 1.0, g, h, d)
        E, V = np.linalg.eigh(off.toarray() + np.diag(diag))
        psi0, phi = ed.initial_states(L, 0); c0 = V.T @ psi0; c1 = V.T @ phi
        y = []; x = []
        for (t0, t1) in WINDOWS:
            acc = np.zeros(2)
            for t in np.random.default_rng(1).uniform(t0, t1, NPT):
                ph = np.exp(-1j * E * t); acc += ed.qfi_site(V @ (ph * c0), V @ (ph * c1), L, 0)
            y.append(acc[0] / NPT); x.append(acc[1] / NPT)
        rows.append([d] + y + x); print(f"L={L} g={g} h={h} d={d:+.3f} Fy " + " ".join(f"{v:.3f}" for v in y), flush=True)
    np.savetxt(f, rows, fmt="%.6e", header=f"L={L} gamma={g} h={h} J=1; full diagonalization; windows tJ: {WINDOWS}; {NPT} random times per window\n"
               "delta | time-averaged F_y in each window | time-averaged F_x in each window")

if __name__ == "__main__":
    for (L, g, h) in ((10, 1.0, 0.2), (12, 1.0, 0.2), (12, 1.0, 0.4), (12, 1.0, 0.0), (10, 1.0, 0.4), (12, 0.5, 0.2)):
        scan(L, g, h)
