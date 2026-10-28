"""
bdg.py -- reference implementation of the free-fermion (delta = 0) XY chain in the full Majorana basis.

Model (Pauli matrices, open boundaries, units J = 1 unless J is passed):
    H = -(J/2) sum_j [(1+gam) sx_j sx_{j+1} + (1-gam) sy_j sy_{j+1}] - sum_j h_j sz_j .
Jordan-Wigner Majoranas  a_j = (prod_{l<j} sz_l) sx_j,  b_j = (prod_{l<j} sz_l) sy_j,  ordered as
    g = (a_1, ..., a_L, b_1, ..., b_L)   (index j for a_j, L + j for b_j, 0-based),
so that H = (i/4) g^T A g with A real antisymmetric (2L x 2L).  The Heisenberg equation gives
    g(t) = G(t) g(0),   G(t) = exp(A t), real orthogonal.

Functions
    Amat(L, J, gam, h)  -> A  (h scalar or array of length L)
    spec(A)             -> (w, v): eigen-decomposition of the Hermitian matrix iA; G(t) = sum_n exp(-i w_n t) v_n v_n^+
    Gt(w, v, t)         -> G(t) (real 2L x 2L)
    plateaus(w, v, L, k, j, wc) -> infinite-time quantities with modes |w| < wc treated as zero modes
                           (not used for the paper data; kept as a reference)

The fast code used for the production runs is free_fermion.py (chiral L x L block); bdg.py serves as the
independent 2L x 2L implementation against which it is tested (tests/test_free.py, F2) and as the reference
propagator for the interacting code at delta = 0 (tests/test_interacting.py, T2).
No input files, no output files.  Imported by free_fermion.py, check_edge_cases.py and the tests.
"""
import numpy as np
def Amat(L, J, gam, h):
    """H=(i/4) g^T A g ; ordering a_1..a_L, b_1..b_L ; h array"""
    h = np.broadcast_to(np.asarray(h, float), (L,))
    A = np.zeros((2*L, 2*L))
    a = lambda j: j; b = lambda j: L + j
    for j in range(L):
        # field: -h_j sz_j = -h_j (-i a_j b_j)  ->  A[a_j, b_j] = 2 h_j
        A[a(j), b(j)] += 2*h[j]; A[b(j), a(j)] -= 2*h[j]
        if j < L-1:
            # exchange: sx_j sx_{j+1} = -i b_j a_{j+1},  sy_j sy_{j+1} = i a_j b_{j+1}
            A[b(j), a(j+1)] += J*(1+gam); A[a(j+1), b(j)] -= J*(1+gam)
            A[a(j), b(j+1)] += -J*(1-gam); A[b(j+1), a(j)] -= -J*(1-gam)
    return A
def spec(A):
    w, v = np.linalg.eigh(1j*A)      # G(t)=sum_n exp(-i w_n t)|n><n|
    return w, v
def Gt(w, v, t):
    return ((v*np.exp(-1j*w*t)) @ v.conj().T).real
def plateaus(w, v, L, k=0, j=0, wc=1e-6):
    """returns dict: avg|W|^2 (modes |w|<wc treated as degenerate/stationary) and DC |Wbar|^2, for Ry and Rx"""
    out = {}
    low = np.abs(w) < wc
    for name, src in (('y', k), ('x', L+k)):
        tot_avg = 0.0; dc = 0.0
        for dst in (j, L+j):
            amp = v[dst, :]*v[src, :].conj()          # P_n(dst,src)
            s0 = amp[low].sum()
            tot_avg += (np.abs(amp[~low])**2).sum() + abs(s0)**2
            dc += abs(s0)**2
        out['avg_'+name] = tot_avg; out['dc_'+name] = dc
    return out
