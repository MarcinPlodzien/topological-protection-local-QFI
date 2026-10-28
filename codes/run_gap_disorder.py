"""Gap of the disordered chain, gamma = 1, h = 0.5J, box fields h_j in [h - Delta, h + Delta].
The single-particle energies are the singular values s_n of B (codes/free_fermion.py); the lowest one is the
edge-mode splitting eps_0 in the topological phase, so the bulk gap is the second smallest.
Output: data/gap/gap_h0.50.txt with median, 10th percentile and minimum over realizations.

Delta = 0, 0.25, 0.4, 0.6, 1, 1.5, 2; L = 100, 200, 400, 800; 100 realizations (one at Delta = 0); one generator
np.random.default_rng(7) drawn sequentially over all (Delta, L).  J = 1, delta = 0.
Output columns: Delta L median p10 min of the second smallest singular value of B (data/gap/gap_h0.50.txt).
Quoted in the reply (Referee A, point 4).
Run (from the repository root):  python3 codes/run_gap_disorder.py      (single process)
Runtime 79 s, peak memory 47 MB (2 threads).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from free_fermion import Bmat

J, gam, h, nreal = 1.0, 1.0, 0.5, 100
rng = np.random.default_rng(7)
rows = []
for D in (0.0, 0.25, 0.4, 0.6, 1.0, 1.5, 2.0):
    for L in (100, 200, 400, 800):
        g = []
        for _ in range(nreal if D > 0 else 1):
            s = np.linalg.svd(Bmat(L, J, gam, h + rng.uniform(-D, D, L)), compute_uv=False)
            g.append(np.sort(s)[1])
        g = np.array(g); rows.append([D, L, np.median(g), np.percentile(g, 10), g.min()])
        print('Delta=%.2f L=%4d  median %.3e  p10 %.3e  min %.3e   bound 2(J-h-Delta)=%.2f' % (D, L, rows[-1][2], rows[-1][3], rows[-1][4], 2 * (J - h - D)), flush=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'gap', 'gap_h0.50.txt')
np.savetxt(out, np.array(rows), header='gamma=1 h=0.5 J=1, %d realizations; Delta L median p10 min of second smallest singular value of B' % nreal, fmt='%.6e')
