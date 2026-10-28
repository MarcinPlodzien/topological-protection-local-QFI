#!/usr/bin/env python3
"""
run_interacting.py -- production runs for the interacting chain (FigR9, FigR10, TableR2).

Writes data/interacting/trace_L{L}_g{gamma}_h{h}_d{delta}.txt with columns
    tJ   F_y(t)   F_x(t)
(single-site QFI of the edge spin, encoding and read-out at site 1, theta = 0).
Existing files are skipped, so the script can be restarted.

    python3 codes/run_interacting.py [nproc]

Model: H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z + delta sum ZZ (Pauli), OBC, J = 1; initial state all
down; R_y and R_x at site 1; QFI of site 1 at theta = 0 (ed_interacting.qfi_site).
Propagator: ed_interacting.cheb_traces_parity (Chebyshev expansion of exp(-iHt) in the two parity sectors, series
cut at |J_k| < 1e-12), output every dt = 0.5/J (nsub = 4 outputs per expansion, i.e. one expansion per 2/J).
The header of every file records the final norm error (at most 1.0e-12 in modulus over all traces).
Task list (64 traces; tmax = 1000/J unless stated):
    gamma = 1, h = 0.2J and 0.4J, delta = 0.05, 0.1, 0.15, 0.2, 0.3, 0.4:  L = 12, 14, 16, 18
    gamma = 1, L = 20, tmax = 300/J:  (h, delta) = (0.2, 0.2), (0.2, 0.3), (0.4, 0.1), (0.4, 0.2)
    EXTRA = (gamma, h, delta) = (0.5, 0.2, 0.1), (0.5, 0.2, 0.2), (1, 0.2, -0.2), (1, 0.4, -0.1):  L = 12, 14, 16
The L = 12 and 14 traces were computed with the full-space propagator ed_interacting.cheb_traces, which gives the
same result (test T7).  The XYZ-line traces (h = 0) come from run_xyz_line.py.
Runtime per trace with 12 processes on a 16-core machine:
L = 16: 1450-2000 s; L = 18: 5600-10200 s; L = 20 (tJ <= 300): 7100-13200 s.  Each process uses one thread.
Peak memory per process: 0.28 GB at L = 18 and 0.93 GB at L = 20 (runs to tJ = 4).
The L = 20 traces in data/interacting/ (to tJ = 1000) are written by run_interacting_long.py; since existing files are
skipped, this script does not overwrite them.
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ed_interacting as ed
OUT = os.path.join(HERE, '..', 'data', 'interacting'); os.makedirs(OUT, exist_ok=True)

J = 1.0
# ordered by priority: the first twelve tasks decide the size dependence of the lifetime
tasks = []
PRI = {0.2: (0.1, 0.15, 0.2, 0.3), 0.4: (0.05, 0.1, 0.15, 0.2)}
tasks += [(18, 1.0, h, d, 1000.0) for h in PRI for d in PRI[h]]
tasks += [(20, 1.0, 0.2, 0.3, 300.0), (20, 1.0, 0.4, 0.2, 300.0)]
tasks += [(16, 1.0, h, d, 1000.0) for h in (0.2, 0.4) for d in (0.05, 0.1, 0.15, 0.2, 0.3, 0.4)]
tasks += [(20, 1.0, 0.2, 0.2, 300.0), (20, 1.0, 0.4, 0.1, 300.0)]
EXTRA = ((0.5, 0.2, 0.1), (0.5, 0.2, 0.2), (1.0, 0.2, -0.2), (1.0, 0.4, -0.1))
tasks += [(16, g, h, d, 1000.0) for g, h, d in EXTRA]
tasks += [(18, 1.0, 0.2, 0.05, 1000.0), (18, 1.0, 0.2, 0.4, 1000.0), (18, 1.0, 0.4, 0.3, 1000.0), (18, 1.0, 0.4, 0.4, 1000.0)]
for L in (14, 12):
    tasks += [(L, 1.0, h, d, 1000.0) for h in (0.2, 0.4) for d in (0.05, 0.1, 0.15, 0.2, 0.3, 0.4)]
    tasks += [(L, g, h, d, 1000.0) for g, h, d in EXTRA]

def fname(L, g, h, d):
    return os.path.join(OUT, f"trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt")

def work(task):
    L, g, h, d, tmax = task
    f = fname(L, g, h, d)
    if os.path.exists(f):
        return f"skip {f}"
    t0 = time.time()
    t, r, nrm = ed.cheb_traces_parity(L, J, g, h, d, tmax, dt=0.5, nsub=4)
    np.savetxt(f, np.column_stack([t, r]), fmt="%.6f %.8e %.8e",
               header=f"H = -(J/2)sum[(1+g)XX+(1-g)YY] - h sum Z + delta sum ZZ (Pauli), OBC, J=1\n"
                      f"L={L} gamma={g} h={h} delta={d}; Chebyshev propagation, final norm-1 = {nrm-1:.1e}\n"
                      f"tJ  F_y  F_x   (QFI of site 1, rotation at site 1, theta=0, initial state all down)")
    return f"done L={L} g={g} h={h} d={d} in {time.time()-t0:.0f}s norm-1={nrm-1:.1e}"

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print(f"{len(tasks)} tasks, {nproc} processes", flush=True)
    with Pool(nproc) as p:
        for msg in p.imap_unordered(work, tasks, chunksize=1):
            print(msg, flush=True)
