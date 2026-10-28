"""Level statistics of the interacting chain: is the Hamiltonian of the lifetime runs integrable or chaotic?

Model (as in all interacting runs, Pauli matrices, open boundaries, J = 1):
    H = -(J/2) sum_j [(1+gamma) X_j X_{j+1} + (1-gamma) Y_j Y_{j+1}] - h sum_j Z_j + delta sum_j Z_j Z_{j+1},  gamma = 1.
For delta = 0 (free fermions) and for h = 0 (free fermions in a rotated frame) the chain is integrable; for h != 0 and
delta != 0 it is expected to be chaotic. This script tests that expectation at the parameters of TableR2.

Method.
1. H conserves the parity of the number of up spins. We take the even sector, which contains the initial state of the
   protocol (all spins down). For L = 14 it has 2^13 = 8192 states.
2. With uniform h, J, delta the open chain is also symmetric under the reflection j -> L+1-j. Levels of different
   symmetry sectors do not repel, so mixing them would mimic Poisson statistics. We therefore split the even sector into
   the reflection-even and reflection-odd blocks (4160 and 4032 states for L = 14) and diagonalize each block fully
   (dense, numpy.linalg.eigvalsh).
3. For each block we compute the ratio of consecutive level spacings, r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1}),
   s_n = E_{n+1} - E_n, and average it over the central half of the spectrum (the band edges have few levels and
   non-generic statistics). The ratio needs no unfolding of the density of states.
   Reference values: Poisson (integrable) <r> = 0.386; GOE (chaotic, real symmetric H) <r> = 0.531.
   The case delta = 0 is included as a check of the method (it must give Poisson).

   The block Hamiltonian B^T H B is formed with sparse matrices (each column of B has at most two entries) and
   converted to a dense matrix only for the diagonalization.

Output: one line per (h, delta): h, delta, dimensions of the two reflection blocks, <r> in each block.
Run (from the folder of the revision), L given as argument (default 14):
    python3 codes/check_level_statistics.py 14 > data/interacting/level_statistics_L14.txt
    python3 codes/check_level_statistics.py 16 > data/interacting/level_statistics_L16.txt
An optional second argument n skips the first n parameter sets (to continue an interrupted run):
    python3 codes/check_level_statistics.py 16 3 > data/interacting/level_statistics_L16_part2.txt
Cost: L = 14, blocks of about 4100 states, about 2 min on 8 threads; L = 16, blocks of about 16400 states,
about 2 GB per block and 10-20 min per parameter set.
"""
import os, sys
import numpy as np
import scipy.sparse as sps
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_interacting as ed

L = int(sys.argv[1]) if len(sys.argv) > 1 else 14
FIRST = int(sys.argv[2]) if len(sys.argv) > 2 else 0
SETS = ((0.2, 0.0), (0.2, 0.05), (0.2, 0.1), (0.2, 0.2), (0.2, 0.3), (0.4, 0.05), (0.4, 0.1), (0.4, 0.2))


def refl(s):
    """Basis state (bit string s, bit i = spin i+1, 1 = up) mirrored by j -> L+1-j."""
    r = 0
    for i in range(L):
        if (s >> i) & 1:
            r |= 1 << (L - 1 - i)
    return r


def ratio(E):
    """Mean spacing ratio <r> over the central half of the sorted spectrum E."""
    E = np.sort(E); n = len(E)
    s = np.diff(E[n // 4: 3 * n // 4])
    s = s[s > 1e-12]                         # exact degeneracies would give r = 0; none occur after symmetry resolution
    return float((np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])).mean())


print(f"# gamma=1 L={L} even parity; columns: h delta dim(R=+1) dim(R=-1) <r>(R=+1) <r>(R=-1); Poisson 0.386, GOE 0.531")
for h, d in SETS[FIRST:]:
    # st: the basis states of the even sector (sorted bit strings); H: the sector Hamiltonian (sparse, real)
    st, H = ed.build_H_sector(L, 1.0, 1.0, h, d, 0)
    M = len(st)
    pos = {s: i for i, s in enumerate(st)}
    R = np.array([pos[refl(s)] for s in st])         # reflection as a permutation of the sector basis
    res = []
    for sgn in (1, -1):
        # orthonormal basis of the block with reflection eigenvalue sgn:
        # a mirror-symmetric state alone (only for sgn = +1), or (|s> + sgn |R s>)/sqrt(2) for each pair s != R s
        rows, cidx, vals = [], [], []; seen = set(); c = 0
        for i in range(M):
            if i in seen:
                continue
            j = R[i]; seen.update((i, j))
            if i == j:
                if sgn == 1:
                    rows.append(i); cidx.append(c); vals.append(1.0); c += 1
            else:
                rows += [i, j]; cidx += [c, c]; vals += [1 / np.sqrt(2), sgn / np.sqrt(2)]; c += 1
        B = sps.csr_matrix((vals, (rows, cidx)), shape=(M, c))     # M x (block dimension), orthonormal columns
        Hb = (B.T @ H @ B).toarray()                                 # block Hamiltonian, dense for diagonalization
        E = np.linalg.eigvalsh(Hb); del Hb                           # all eigenvalues of the block
        res.append((len(E), ratio(E)))
    print(f"{h:.2f} {d:.2f} {res[0][0]} {res[1][0]} {res[0][1]:.3f} {res[1][1]:.3f}", flush=True)
