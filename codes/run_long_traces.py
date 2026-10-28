#!/usr/bin/env python3
"""
run_long_traces.py -- edge-QFI traces to tJ = 10^4 for the slowly decaying sets, so that the decay is resolved
within the trace and the lifetime is not extrapolated.

Model, initial state, encoding and read-out as in run_interacting.py (gamma = 1, site 1, theta = 0, all spins down).
Propagator: ed_interacting.cheb_traces_parity_fast (equal to cheb_traces_parity to 1e-14, test T9).
Output: data/interacting_long/trace_L{L}_g1.00_h{h}_d{delta}_t10000.txt, columns tJ, F_y, F_x.
Incremental saving: every 50/J the state vectors and the trace so far are written to <output>.ckpt.npz
(through a temporary file and an atomic rename) and the trace so far to <output>.ckpt.partial.txt.
A restarted run continues from the checkpoint; finished outputs are skipped.

    python3 codes/run_long_traces.py [nproc] [threads per process]

Cost on this machine: about 2.4 s per unit of tJ at L = 20 with two such runs side by side (6 threads each),
i.e. about 7 h per L = 20 trace; L = 18 about a quarter of that, L = 16 about a sixteenth.
Measured with 2 processes x 6 threads: L = 20 22373 s and 22716 s, L = 18 4066 s and 4073 s, L = 16 743 s and 649 s.  Final norm errors 7.5e-14 to 6.8e-12 (file headers).
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import get_context
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, '..', 'data', 'interacting_long'); os.makedirs(OUT, exist_ok=True)

J, GAMMA, TMAX = 1.0, 1.0, 10000.0
# largest sizes first, so that the two long L = 20 runs start together
TASKS = [(20, 0.2, 0.2), (20, 0.4, 0.1), (18, 0.2, 0.2), (18, 0.4, 0.1), (16, 0.2, 0.2), (16, 0.4, 0.1)]

def fname(L, h, d):
    return os.path.join(OUT, f"trace_L{L}_g{GAMMA:.2f}_h{h:.2f}_d{d:.2f}_t{TMAX:.0f}.txt")

def work(task):
    L, h, d, threads = task
    import ed_interacting as ed
    f = fname(L, h, d)
    if os.path.exists(f):
        return f"skip {f}"
    t0 = time.time()
    t, r, nrm = ed.cheb_traces_parity_fast(L, J, GAMMA, h, d, TMAX, dt=0.5, nsub=4, threads=threads,
                                           checkpoint=f[:-4] + '.ckpt.npz', every=25)
    np.savetxt(f + '.tmp', np.column_stack([t, r]), fmt="%.6f %.8e %.8e",
               header=f"H = -(J/2)sum[(1+g)XX+(1-g)YY] - h sum Z + delta sum ZZ (Pauli), OBC, J=1\n"
                      f"L={L} gamma={GAMMA} h={h} delta={d}; Chebyshev propagation (cheb_traces_parity_fast, "
                      f"{threads} threads), final norm-1 = {nrm-1:.1e}\n"
                      f"tJ  F_y  F_x   (QFI of site 1, rotation at site 1, theta=0, initial state all down)")
    os.replace(f + '.tmp', f)
    return f"done L={L} h={h} d={d} to tJ={TMAX:.0f} in {time.time()-t0:.0f}s norm-1={nrm-1:.1e}"

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    tasks = [(L, h, d, threads) for L, h, d in TASKS]
    print(f"{time.strftime('%Y-%m-%d %H:%M')}  {len(tasks)} traces to tJ={TMAX:.0f}, {nproc} processes x {threads} threads",
          flush=True)
    with get_context("spawn").Pool(nproc) as p:
        for msg in p.imap_unordered(work, tasks, chunksize=1):
            print(time.strftime('%Y-%m-%d %H:%M'), msg, flush=True)
    print(time.strftime('%Y-%m-%d %H:%M'), "all done", flush=True)
