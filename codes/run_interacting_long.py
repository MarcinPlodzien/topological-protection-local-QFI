#!/usr/bin/env python3
"""
run_interacting_long.py -- edge-QFI traces of the interacting chain to tJ = 1000 for the largest sizes.

Same model, initial state, encoding and read-out as run_interacting.py (gamma = 1, rotation and read-out at
site 1, theta = 0, all spins down), propagated with ed_interacting.cheb_traces_parity_fast (numba, parallel over
rows; identical to cheb_traces_parity to 1e-14, test T9). Output in data/interacting/:
    trace_L{L}_g{gamma}_h{h}_d{delta}.txt     columns tJ, F_y, F_x
Checkpoints (state vectors and partial trace) are written every 50/J and the run resumes from them;
finished files are skipped.

    python3 codes/run_interacting_long.py L [nproc] [threads per process] [number of sets]
    e.g. python3 codes/run_interacting_long.py 20 2 6       (all twelve sets)
         python3 codes/run_interacting_long.py 22 2 6 8     (the first eight sets)

Cost on this machine (16 cores): L = 20, 6 threads, about 1.2 s per unit of tJ, i.e. 20 min per trace.
Measured with 2 processes x 6 threads: L = 20 2366-7386 s per trace; L = 22 9498-11368 s per trace.
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import get_context
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, '..', 'data', 'interacting'); os.makedirs(OUT, exist_ok=True)   # read by make_figures.py (TableR2, Figs. 4, S4)

J, GAMMA, TMAX = 1.0, 1.0, 1000.0
# the gamma = 1 parameter sets of TableR2, in order of priority (sets whose lifetime still depends on L first)
SETS = [(0.4, 0.1), (0.4, 0.2), (0.2, 0.2), (0.2, 0.3), (0.4, 0.05), (0.4, 0.15), (0.4, 0.3), (0.2, 0.15),
        (0.2, 0.1), (0.2, 0.05), (0.2, 0.4), (0.4, 0.4)]

def fname(L, h, d):
    return os.path.join(OUT, f"trace_L{L}_g{GAMMA:.2f}_h{h:.2f}_d{d:.2f}.txt")

def work(task):
    L, h, d, threads = task
    import ed_interacting as ed
    f = fname(L, h, d)
    if os.path.exists(f):
        return f"skip {f}"
    t0 = time.time()
    t, r, nrm = ed.cheb_traces_parity_fast(L, J, GAMMA, h, d, TMAX, dt=0.5, nsub=4, threads=threads,
                                           checkpoint=f[:-4] + '.ckpt.npz', every=25)
    np.savetxt(f, np.column_stack([t, r]), fmt="%.6f %.8e %.8e",
               header=f"H = -(J/2)sum[(1+g)XX+(1-g)YY] - h sum Z + delta sum ZZ (Pauli), OBC, J=1\n"
                      f"L={L} gamma={GAMMA} h={h} delta={d}; Chebyshev propagation (cheb_traces_parity_fast, "
                      f"{threads} threads), final norm-1 = {nrm-1:.1e}\n"
                      f"tJ  F_y  F_x   (QFI of site 1, rotation at site 1, theta=0, initial state all down)")
    return f"done L={L} h={h} d={d} in {time.time()-t0:.0f}s norm-1={nrm-1:.1e}"

if __name__ == "__main__":
    L = int(sys.argv[1]); nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    nsets = int(sys.argv[4]) if len(sys.argv) > 4 else len(SETS)      # the first nsets sets of SETS
    tasks = [(L, h, d, threads) for h, d in SETS[:nsets]]
    print(f"{time.strftime('%Y-%m-%d %H:%M')}  L={L}: {len(tasks)} traces to tJ={TMAX:.0f}, "
          f"{nproc} processes x {threads} threads", flush=True)
    with get_context("spawn").Pool(nproc) as p:
        for msg in p.imap_unordered(work, tasks, chunksize=1):
            print(time.strftime('%H:%M'), msg, flush=True)
    print(time.strftime('%Y-%m-%d %H:%M'), "all done", flush=True)
