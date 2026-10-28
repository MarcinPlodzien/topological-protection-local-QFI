"""Infinite-temperature autocorrelator of the edge spin, A(t) = 2^-L Tr[sx_1(t) sx_1], by full diagonalization,
for comparison with the edge QFI of the polarized initial state (which equals A^2 only for delta = 0).
Output: data/interacting/autocorr_L{L}.txt.   Usage: python3 codes/run_autocorrelator.py [L=12]

Model: as in run_interacting.py (Pauli, OBC, J = 1, gamma = 1).  Parameter sets (gamma, h, delta) =
(1, 0.2, 0.2), (1, 0.2, -0.2), (1, 0.2, 0.3), (1, 0.4, 0.2).
A_inf(t) = 2^-L sum_{m,n} |<m|sx_1|n>|^2 cos((E_m - E_n) t) from the full spectrum; sx_1 flips bit 0 of the basis index.
A_inf and the edge QFI F_y of the all-down state are averaged over the window [t-5, t] for t = 10, 20, 50, 100, 200,
400, 700, 1000 (A_inf on 11 equally spaced times, F_y on the stored trace points of
data/interacting/trace_L{L}_g.._h.._d...txt, which must exist).  The lifetime of A_inf^2 printed to stdout is
-1/slope of a linear fit of ln A_inf^2 over the eight times (quoted in the reply, and in the reply).
Output columns: gamma h delta tJ A_inf A_inf^2 F_y.
Runtime 189 s for L = 12, peak memory 1.1 GB (2 threads).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_interacting as ed

L = int(sys.argv[1]) if len(sys.argv) > 1 else 12
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'interacting')
TS = np.array([10, 20, 50, 100, 200, 400, 700, 1000.])
PAR = [(1.0, 0.2, 0.2), (1.0, 0.2, -0.2), (1.0, 0.2, 0.3), (1.0, 0.4, 0.2)]
rows = []
for g, h, d in PAR:
    N = 1 << L
    diag, off = ed.build_H(L, 1.0, g, h, d)
    E, V = np.linalg.eigh(off.toarray() + np.diag(diag))
    M = (V.T @ V[np.arange(N) ^ 1, :]) ** 2            # |<m|sx_1|n>|^2
    dE = E[:, None] - E[None, :]
    A = np.array([np.mean([(M * np.cos(dE * x)).sum() / N for x in np.linspace(t - 5, t, 11)]) for t in TS])
    tr = np.loadtxt(os.path.join(D, f'trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt'))
    F = np.array([tr[(tr[:, 0] >= t - 5 - 1e-9) & (tr[:, 0] <= t + 1e-9), 1].mean() for t in TS])
    tauA2 = -1.0 / np.polyfit(TS, np.log(A ** 2), 1)[0]
    print(f'gamma={g} h={h} delta={d:+.2f}: lifetime of A^2 = {tauA2:.3g}', flush=True)
    for t, a, f in zip(TS, A, F):
        rows.append((g, h, d, t, a, a * a, f))
np.savetxt(os.path.join(D, f'autocorr_L{L}.txt'), rows, fmt='%.2f %.2f %+.2f %7.1f %.6f %.6f %.6f',
           header=f'L={L}; full diagonalization; averages over the window [t-5, t]\n'
                  'gamma h delta tJ A_inf A_inf^2 F_y(edge QFI, all-down initial state)')
