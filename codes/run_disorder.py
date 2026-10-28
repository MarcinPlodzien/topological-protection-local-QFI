#!/usr/bin/env python3
"""
run_disorder.py -- edge QFI in the chain with random fields h_j = h + uniform[-Delta, Delta].

For every realization: stationary-state QFI and time-averaged QFI (window [0, T]) for R_y and R_x,
encoding and read-out at site 1, plus the zero-mode Lyapunov exponent of that realization.
Writes data/disorder/disorder_L{L}_g{gamma}_h{h}_T{T}.txt with one row per Delta:
    Delta  lyap_mean | stat_y: mean median q25 q75 | stat_x: ... | avg_y: ... | avg_x: ... | frac(stat_y > 10 stat_x)
The chain is topological iff the disorder-averaged Lyapunov exponent is negative
(gamma = 1:  <ln|h_j|> < ln J).

    python3 codes/run_disorder.py [nproc]

Model: free XY chain (delta = 0), J = 1, k = j = 1, window [0, T], fields h_j = h + uniform[-Delta, Delta],
Delta = 0, 0.25, ..., 4 (17 values; at Delta = 0 one realization).  Runs (L, gamma, h, T, realizations), see RUNS:
    (200,1,0.5,1e3,1000) (100,1,0.5,1e3,1000) (400,1,0.5,1e3,300) (200,1,0.5,1e2,1000) (200,1,0.5,1e4,1000)
    (200,0.5,0.5,1e3,1000) (200,1,1.5,1e3,1000) (200,1,1.5,1e2,1000) (200,1,1.5,1e4,1000)
Seeds: realization s at the iD-th Delta uses np.random.default_rng(1000003*iD + s + 17*L + int(1000*gamma) + int(100*h)).
Output columns (data/disorder/disorder_L{L}_g{gamma}_h{h}_T{T}.txt, %.6e):
    0 Delta | 1 mean Lyapunov exponent | 2-5 stat_y mean, median, q25, q75 | 6-9 stat_x | 10-13 avg_y | 14-17 avg_x
    | 18 fraction of realizations with stat_y > 10 stat_x
Also, before the finite chains, the exact semi-infinite reference (gamma = 1, T -> infinity, 2e4 sites, 4000 realizations,
seed 7000003*iD + s + int(1000*h)): data/disorder/semi_infinite_g1.00_h1.50.txt,
columns Delta | mean median q25 q75 of I_st = [sum_j prod_{l<j}(h_l/J)^2]^-2 | fraction with I_st > 1e-6;
Delta grid 1..3 step 0.05.  The reference for h = 0.5 (40000 realizations) is written by run_semi_infinite_refine.py.
Existing files are skipped, so the script can be restarted and extended.
Runtime: 1-2 h for all tables with 4 processes; the two h = 1.5J tables with TJ = 1e2 and 1e4 take
about 3 min with 4 processes.  Memory: below 1 GB per process.
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from free_fermion import chiral_modes, late_time, lyapunov_zero_mode, stationary_exact_ising
OUT = os.path.join(HERE, '..', 'data', 'disorder'); os.makedirs(OUT, exist_ok=True)
DELTAS = np.arange(0, 4.01, 0.25)
#        L   gamma  h    T      realizations
RUNS = [(200, 1.0, 0.5, 1e3, 1000), (100, 1.0, 0.5, 1e3, 1000), (400, 1.0, 0.5, 1e3, 300),
        (200, 1.0, 0.5, 1e2, 1000), (200, 1.0, 0.5, 1e4, 1000),
        (200, 0.5, 0.5, 1e3, 1000), (200, 1.0, 1.5, 1e3, 1000),
        (200, 1.0, 1.5, 1e2, 1000), (200, 1.0, 1.5, 1e4, 1000)]      # dependence on T for the chain that is trivial without disorder

def one(args):
    L, g, h, T, D, seed = args
    rng = np.random.default_rng(seed + 17 * L + int(1000 * g) + int(100 * h))
    hj = h + rng.uniform(-D, D, L)
    r = late_time(*chiral_modes(L, 1.0, g, hj), T)
    return [r['stat_y'], r['stat_x'], r['avg_y'], r['avg_x'], lyapunov_zero_mode(hj, 1.0, g)]

def semi_infinite(args):
    h, D, seed = args
    return stationary_exact_ising(h + np.random.default_rng(seed).uniform(-D, D, 20000))

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(nproc) as p:
        # exact reference for gamma = 1: semi-infinite chain, T -> infinity (zero-mode weight on the first site)
        for hs in (1.5,):                     # h = 0.5: codes/run_semi_infinite_refine.py (40000 realizations)
            f = os.path.join(OUT, f'semi_infinite_g1.00_h{hs:.2f}.txt')
            if os.path.exists(f): continue
            Dg = np.arange(0, 4.01, 0.125) if hs == 0.5 else np.arange(1.0, 3.001, 0.05); rows = []
            for iD, D in enumerate(Dg):
                a = np.array(p.map(semi_infinite, [(hs, D, 7000003 * iD + s + int(1000 * hs)) for s in range(4000 if D > 0 else 1)], chunksize=50))
                rows.append([D, a.mean(), np.median(a), np.quantile(a, 0.25), np.quantile(a, 0.75), np.mean(a > 1e-6)])
            np.savetxt(f, rows, fmt='%.6e', header=f'gamma=1 h={hs}, semi-infinite chain (2e4 sites), T->infinity, 4000 realizations\nDelta | mean median q25 q75 of I_st=[sum_j prod_{{l<j}}(h_l/J)^2]^-2 | fraction of realizations with I_st>1e-6')
            print('wrote', f, flush=True)
        for (L, g, h, T, nr) in RUNS:
            f = os.path.join(OUT, f"disorder_L{L}_g{g:.2f}_h{h:.2f}_T{T:.0e}.txt")
            if os.path.exists(f): continue
            rows = []
            for iD, D in enumerate(DELTAS):
                n = nr if D > 0 else 1
                a = np.array(p.map(one, [(L, g, h, T, D, 1000003 * iD + s) for s in range(n)], chunksize=8))
                row = [D, a[:, 4].mean()]
                for c in range(4):
                    row += [a[:, c].mean(), np.median(a[:, c]), np.quantile(a[:, c], 0.25), np.quantile(a[:, c], 0.75)]
                row.append(np.mean(a[:, 0] > 10 * a[:, 1]))
                rows.append(row); print(f"L={L} g={g} h={h} T={T:.0e} " + " ".join(f"{x:.4f}" for x in row), flush=True)
            np.savetxt(f, rows, fmt="%.6e", header=f"L={L} gamma={g} h={h} J=1, window [0,T] with T={T}, {nr} realizations, h_j=h+U[-Delta,Delta]\n"
                       "Delta lyap | stat_y mean med q25 q75 | stat_x mean med q25 q75 | avg_y mean med q25 q75 | avg_x mean med q25 q75 | frac(stat_y>10 stat_x)")
