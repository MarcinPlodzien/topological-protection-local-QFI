"""Edge QFI (R_y, R_x) of the interacting chain with a trivial edge state:
gamma = 1, h_j = 1.8 J for j >= 2, h_1 = 0 (zero-energy edge state) or 0.3 J (finite-energy bound state),
delta * sum Z_j Z_{j+1}, all spins down, rotation and read-out of the edge spin. Exact propagation in the parity
sectors (ed_interacting.cheb_traces_parity), output every 0.5/J to tJ = 1000."""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "codes"))
import ed_interacting as ed
OUT = os.path.dirname(os.path.abspath(__file__))
def run(task):
    L, h1, d = task
    f = os.path.join(OUT, f"trivial_L{L}_h1{h1:.2f}_d{d:.2f}.txt")
    if os.path.exists(f): return f
    h = np.full(L, 1.8); h[0] = h1
    t0 = time.time()
    t, F, nrm = ed.cheb_traces_parity(L, 1.0, 1.0, h, d, 1000.0, dt=0.5, nsub=4, k=0, j=0)
    np.savetxt(f, np.column_stack([t, F]), fmt="%.6f %.9e %.9e",
               header=f"gamma=1 h_j=1.8 (j>=2) h_1={h1} delta={d} L={L}; norm={nrm:.12f}; runtime {time.time()-t0:.0f} s\ntJ F_y F_x")
    return f
if __name__ == "__main__":
    tasks = [(L, h1, d) for L in (14, 16) for (h1, d) in [(0.0, 0.0), (0.0, 0.05), (0.0, 0.1), (0.0, 0.2), (0.3, 0.1), (0.3, 0.2)]]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        for f in p.imap_unordered(run, tasks): print("done", f, flush=True)
