"""Interacting integrable line h = 0, gamma != 1 (XYZ chain, exact strong zero mode): edge QFI by exact propagation.
The plateau must not decay in the thermodynamic limit; at finite L it decays through the splitting of the two edge modes.
Output: data/interacting/trace_L{L}_g0.50_h0.00_d0.20.txt (same format as run_interacting.py).
Usage: python3 codes/run_xyz_line.py

Model: H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] + delta sum ZZ with gamma = 0.5, h = 0, delta = 0.2J (Pauli, OBC,
J = 1); this is the XYZ chain, integrable with an exact strong zero mode.  Initial state all down, k = j = 1.
L = 12, 14, 16, 18 to tJ = 1000 and L = 20 to tJ = 300, output every 0.5/J (cheb_traces_parity).
Checkpoints: trace_..ckpt.npz and trace_..ckpt.partial.txt every 50/J; a restarted run continues from them and
removes them at the end.  Existing traces are skipped.
Output columns: tJ F_y F_x (same header as run_interacting.py).  Numbers quoted in the reply (Referee B, point 3); not plotted.
Runtime: L = 18 3578 s and L = 20 3563 s with OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = 4.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_interacting as ed

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'interacting')
G, H, DL = 0.5, 0.0, 0.2
for L, tmax in ((12, 1000), (14, 1000), (16, 1000), (18, 1000), (20, 300)):
    f = os.path.join(D, f'trace_L{L}_g{G:.2f}_h{H:.2f}_d{DL:.2f}.txt')
    if os.path.exists(f):
        continue
    t0 = time.time()
    t, r, nrm = ed.cheb_traces_parity(L, 1.0, G, H, DL, tmax, dt=0.5, nsub=4, progress=25,
                                      checkpoint=f[:-4] + '.ckpt.npz', every=25)   # saved every 50/J
    np.savetxt(f, np.column_stack([t, r]), fmt="%.6f %.8e %.8e",
               header=f"H = -(J/2)sum[(1+g)XX+(1-g)YY] - h sum Z + delta sum ZZ (Pauli), OBC, J=1\n"
                      f"L={L} gamma={G} h={H} delta={DL}; Chebyshev propagation, final norm-1 = {nrm - 1:.1e}\n"
                      f"tJ  F_y  F_x   (QFI of site 1, rotation at site 1, theta=0, initial state all down)")
    print(f'L={L} done in {time.time() - t0:.0f} s', flush=True)
