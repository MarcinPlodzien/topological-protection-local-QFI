#!/usr/bin/env python3
"""run_phase_diagram.py -- stationary-state and time-averaged edge QFI over the (h/J, gamma) plane (FigR1a).
L = 200, k = j = 1, window [0, T] with T = 200/J.  Writes data/free/phase_diagram_*.txt

Model: free XY chain (delta = 0), H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z, OBC, J = 1,
initial state all down, theta = 0; free_fermion.late_time gives the window averages in closed form.
Grid: h/J = linspace(0, 2, 81) (columns), gamma = linspace(0, 1, 41) (rows).
Outputs (data/free/):
  phase_diagram_h.txt, phase_diagram_gamma.txt   the two grids
  phase_diagram_{stat_y,stat_x,avg_y,avg_x}.txt  41 x 81 tables (rows gamma, columns h); stat_y is FigR1(a)
  phase_diagram_exact.txt                        stationary QFI of the clean semi-infinite chain (plateau_exact)
Prints the largest |stat_y - exact| for h < 0.9J, gamma >= 0.2.
Run (from the repository root):  python3 codes/run_phase_diagram.py      (Pool of 2 processes, 1 thread each)
Deterministic.  Runtime 32 s with the pool of 2 processes.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from free_fermion import chiral_modes, late_time, plateau_exact
OUT = os.path.join(HERE, '..', 'data', 'free')
hs = np.linspace(0.0, 2.0, 81); gs = np.linspace(0.0, 1.0, 41); L = 200; T = 200.0
def row(g):
    out = []
    for h in hs:
        r = late_time(*chiral_modes(L, 1.0, g, h), T)
        out.append([r['stat_y'], r['stat_x'], r['avg_y'], r['avg_x'], plateau_exact(h, 1.0, g)])
    return out
if __name__ == "__main__":
    with Pool(2) as p: res = np.array(p.map(row, gs))
    hd = f"rows: gamma grid (phase_diagram_gamma.txt), columns: h/J grid (phase_diagram_h.txt). L={L}, k=j=1, window [0,T], T={T}\n"
    np.savetxt(os.path.join(OUT, "phase_diagram_h.txt"), hs); np.savetxt(os.path.join(OUT, "phase_diagram_gamma.txt"), gs)
    for i, nm in enumerate(("stat_y", "stat_x", "avg_y", "avg_x", "exact")):
        np.savetxt(os.path.join(OUT, f"phase_diagram_{nm}.txt"), res[:, :, i], fmt="%.6e", header=hd + nm)
    m = hs < 0.9
    print("max |stat_y - exact| for h<0.9J, gamma>=0.2:", np.abs(res[8:, m, 0] - res[8:, m, 4]).max())
