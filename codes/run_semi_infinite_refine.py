"""Exact mean stationary edge QFI of the semi-infinite disordered chain (gamma = 1, h = 0.5J, T -> infinity),
I_st = [sum_j prod_{l<j}(h_l/J)^2]^-2, with 40000 realizations per Delta.
Same Delta grid and chain length (2e4 sites) as codes/run_disorder.py. Writes
data/disorder/semi_infinite_g1.00_h0.50.txt.

Model: transverse-field Ising chain (gamma = 1, delta = 0), fields h_j = 0.5 + uniform[-Delta, Delta], J = 1.
For gamma = 1 the left zero mode obeys phi_{j+1} = (h_j/J) phi_j, so the stationary edge QFI (R_y, k = j = 1) of the
semi-infinite chain is I_st = [sum_j prod_{l<j} (h_l/J)^2]^-2 (zero if the sum diverges); the sum is evaluated in
log space (log-sum-exp) to avoid overflow.  Delta = 0, 0.125, ..., 4 (33 values), 40000 realizations per Delta,
drawn in chunks of 400 from a single generator np.random.default_rng(20261014).
Output columns: Delta | mean median q25 q75 of I_st | fraction of realizations with I_st > 1e-6.  Used in FigR6(a), FigR7(a,b).
Run (from the repository root):  python3 codes/run_semi_infinite_refine.py
Runtime: below 30 min with two threads.
Memory: about 0.4 GB (two 400 x 20000 float arrays per chunk plus temporaries).
"""
import os
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'disorder')
f = os.path.join(OUT, 'semi_infinite_g1.00_h0.50.txt')
h, J, NS, N, CH = 0.5, 1.0, 20000, 40000, 400
rng = np.random.default_rng(20261014); rows = []
for D in np.arange(0, 4.01, 0.125):
    if D == 0:
        a = np.array([(1 - h ** 2) ** 2])
    else:
        a = []
        for _ in range(N // CH):
            lr = 2.0 * np.log(np.abs(h + rng.uniform(-D, D, (CH, NS))) / J + 1e-300)
            c = np.concatenate((np.zeros((CH, 1)), np.cumsum(lr, axis=1)[:, :-1]), axis=1)
            m = c.max(axis=1, keepdims=True)
            a.append(np.exp(-2.0 * (m[:, 0] + np.log(np.exp(c - m).sum(axis=1)))))
        a = np.concatenate(a)
    rows.append([D, a.mean(), np.median(a), np.quantile(a, 0.25), np.quantile(a, 0.75), np.mean(a > 1e-6)])
    print('Delta=%.3f mean %.5f +- %.5f median %.4g' % (D, a.mean(), a.std() / np.sqrt(len(a)), np.median(a)), flush=True)
np.savetxt(f, rows, fmt='%.6e', header='gamma=1 h=0.5, semi-infinite chain (2e4 sites), T->infinity, %d realizations\nDelta | mean median q25 q75 of I_st=[sum_j prod_{l<j}(h_l/J)^2]^-2 | fraction of realizations with I_st>1e-6' % N)
print('wrote', f)
