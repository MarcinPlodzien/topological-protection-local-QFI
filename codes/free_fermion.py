"""
free_fermion.py -- exact single-site QFI in the (inhomogeneous) XY chain, delta = 0.

Majorana form  H = (i/4) g^T A g,  g = (a_1..a_L, b_1..b_L),
    a_j = (prod_{l<j} sz_l) sx_j ,   b_j = (prod_{l<j} sz_l) sy_j ,
so that g(t) = G(t) g(0) with G(t) = exp(A t) real orthogonal (see bdg.py).

Rotation exp(-i theta s^a_k/2) of the state |dd...d>, read-out of spin j, theta = 0:
    F_y(t) = G[a_j,a_k]^2 + G[b_j,a_k]^2 ,     F_x(t) = G[a_j,b_k]^2 + G[b_j,b_k]^2 .
(these equal |U+V|^2 and |U-V|^2 of the reply).

Two late-time quantities are defined for the window t in [0, T]:
  * time-averaged QFI            <F>      = mean_t F(t)
  * QFI of the time-averaged state  F[rho_bar]: at theta=0 the Bloch vector of spin j is along z
    and its theta-derivative lies in the xy plane, so F[rho_bar] = (mean_t G[a_j,.])^2 + (mean_t G[b_j,.])^2.
    Only modes with |omega| < ~1/T contribute: it is the zero-frequency (stationary) part of the signal.

Model and conventions.  H = -(J/2) sum_j [(1+gam) sx_j sx_{j+1} + (1-gam) sy_j sy_{j+1}] - sum_j h_j sz_j,
Pauli matrices, open boundaries, J = 1 in all production runs; h may be a scalar or an array of length L
(inhomogeneous or random fields).  Initial state all spins down; rotation exp(-i theta s^a_k/2), a = y (R_y) or
x (R_x); all QFI values at theta = 0.  Sites are 0-based in the code (k = 0 is site 1 of the paper).
In the ordering (a, b) the matrix A has the chiral form [[0, B], [-B^T, 0]]; the singular value decomposition
B = U diag(s) V^T gives the single-particle energies s_n >= 0 and all blocks of G(t) (chiral_modes).
The edge zero mode of the topological phase is the column of U (a-type, left edge) with s_n -> 0.

Main functions
    window_kernel(x, T)            (1/T) int_0^T exp(-i x t) dt, closed form used for all time averages
    Bmat, chiral_modes             chiral block B and its SVD
    late_time(U, s, V, T, k, j)    stationary and time-averaged QFI for the window [0, T] (fast path, used in production)
    late_time_general(w, v, ...)   the same from the 2L x 2L spectrum of bdg.py (reference, slow)
    qfi_t(w, v, L, t, k, j)        F_y(t), F_x(t) without averaging (edge read-out formula)
    lyapunov_zero_mode(h, J, gam)  inverse localization length of the left zero mode (transfer matrices, QR)
    plateau_exact(h, J, gam)       stationary edge QFI of the clean semi-infinite chain
    pfaffian(A), spin_qfi_profile  exact single-SPIN QFI at any read-out site (Jordan-Wigner string included)
    stationary_exact_ising(h, J)   closed form of the stationary edge QFI of a disordered gamma = 1 chain
No input or output files; imported by the run_*.py and check_*.py scripts and by the tests.
"""
import numpy as np
from bdg import Amat, spec, Gt


def window_kernel(x, T):
    """(1/T) int_0^T exp(-i x t) dt"""
    return np.exp(-0.5j * x * T) * np.sinc(x * T / (2 * np.pi))


def Bmat(L, J, gam, h):
    """Chiral block: A = [[0, B], [-B^T, 0]] in the ordering (a_1..a_L, b_1..b_L)."""
    h = np.broadcast_to(np.asarray(h, float), (L,))
    B = np.diag(2.0 * h)
    i = np.arange(L - 1)
    B[i + 1, i] = -J * (1 + gam)
    B[i, i + 1] = -J * (1 - gam)
    return B


def chiral_modes(L, J, gam, h):
    """B = U diag(s) V^T.  Then G_aa = U cos(st) U^T, G_bb = V cos(st) V^T,
       G_ab = U sin(st) V^T, G_ba = -V sin(st) U^T; the mode frequencies are +-s_n."""
    U, s, Vt = np.linalg.svd(Bmat(L, J, gam, h))
    return U, s, Vt.T


def late_time(U, s, V, T, k=0, j=0):
    """Window [0, T].  Returns dict with
         stat_y, stat_x : QFI of the time-averaged state (stationary-state QFI)
         avg_y,  avg_x  : time-averaged QFI
       for a rotation about y or x at site k and read-out of the edge spin j = 0 (for j > 0 the
       expressions give the fermionic-mode quantity; see spin_qfi_profile for the spin QFI)."""
    # cbar_n = mean_t cos(s_n t), sbar_n = mean_t sin(s_n t) over [0, T]
    g = window_kernel(s, T); cbar, sbar = g.real, -g.imag
    # mean_t cos(s_n t) cos(s_m t) = (Cm + Cp)/2 and mean_t sin sin = (Cm - Cp)/2, with Cm, Cp the kernels at s_n -+ s_m
    Cm = window_kernel(s[:, None] - s[None, :], T).real
    Cp = window_kernel(s[:, None] + s[None, :], T).real
    def sq_cos(A): return 0.5 * (A @ (Cm + Cp) @ A)
    def sq_sin(A): return 0.5 * (A @ (Cm - Cp) @ A)
    Aaa = U[j] * U[k]; Aba = V[j] * U[k]; Aab = U[j] * V[k]; Abb = V[j] * V[k]
    return {'stat_y': (Aaa @ cbar) ** 2 + (Aba @ sbar) ** 2, 'avg_y': sq_cos(Aaa) + sq_sin(Aba),
            'stat_x': (Aab @ sbar) ** 2 + (Abb @ cbar) ** 2, 'avg_x': sq_sin(Aab) + sq_cos(Abb)}


def late_time_general(w, v, L, T, k=0, j=0):
    """Same as late_time, from the full 2L x 2L spectrum (slow reference implementation)."""
    out = {}
    g = window_kernel(w, T)
    K = window_kernel(w[:, None] - w[None, :], T)
    for name, src in (('y', k), ('x', L + k)):
        stat = 0.0; avg = 0.0
        for dst in (j, L + j):
            amp = v[dst, :] * v[src, :].conj()
            stat += abs((amp * g).sum()) ** 2
            avg += np.real(amp @ K @ amp.conj())
        out['stat_' + name] = stat; out['avg_' + name] = avg
    return out


def qfi_t(w, v, L, t, k=0, j=0):
    G = Gt(w, v, t)
    return G[j, k] ** 2 + G[L + j, k] ** 2, G[j, L + k] ** 2 + G[L + j, L + k] ** 2


def lyapunov_zero_mode(h, J, gam, rng=None):
    """Largest Lyapunov exponent of the zero-energy recursion for the left (a-type) Majorana,
       J(1+gam) phi_{j+1} = 2 h_j phi_j - J(1-gam) phi_{j-1}   (row b_j of A phi = 0).
       The left end hosts a normalizable zero mode (topological) iff the exponent is negative.
       For gam = 1 it reduces to <ln|h_j|> - ln J."""
    h = np.asarray(h, float)
    if gam == 1.0:
        return np.mean(np.log(np.abs(h) + 1e-300)) - np.log(J)
    Q = np.eye(2); lam = 0.0
    for hj in h:
        M = np.array([[2 * hj / (J * (1 + gam)), -(1 - gam) / (1 + gam)], [1.0, 0.0]])
        Q, R = np.linalg.qr(M @ Q)
        lam += np.log(abs(R[0, 0]))
    return lam / len(h)


def plateau_exact(h, J, gam):
    """Stationary edge QFI (R_y, k=j=1) of the clean semi-infinite chain: |phi_1|^4 for the
       normalized left Majorana mode phi_j ~ (r_+^j - r_-^j)/(r_+ - r_-), (1+gam) r^2 - 2(h/J) r + (1-gam) = 0.
       Zero in the trivial phase.  gam = 1:  (1 - h^2/J^2)^2."""
    x = h / J
    if abs(x) >= 1 or gam <= 0:
        return 0.0
    rp, rm = np.roots([1 + gam, -2 * x, 1 - gam]).astype(complex) if gam != 1 else (complex(x), 0j)
    j = np.arange(1, 4000)
    phi = ((rp ** j - rm ** j) / (rp - rm)).real if abs(rp - rm) > 1e-9 else (j * rp ** (j - 1)).real
    return float((phi[0] ** 2 / (phi ** 2).sum()) ** 2)


# ----------------------------------------------------------------------------------------------
#  Single-SPIN QFI at an arbitrary read-out site j (Jordan-Wigner string included)
# ----------------------------------------------------------------------------------------------
def pfaffian(A):
    """Pfaffian of a complex antisymmetric matrix (Parlett-Reid with pivoting)."""
    A = np.array(A, dtype=complex); n = A.shape[0]
    if n % 2: return 0.0
    # Each step eliminates rows/columns k, k+1 by a congruence transformation that keeps the matrix antisymmetric;
    # Pf(A) is the product of the pivots A[k, k+1], with a sign change for every row/column exchange.
    pf = 1.0 + 0j
    for k in range(0, n - 1, 2):
        kp = k + 1 + np.abs(A[k + 1:, k]).argmax()
        if kp != k + 1:
            A[[k + 1, kp], :] = A[[kp, k + 1], :]; A[:, [k + 1, kp]] = A[:, [kp, k + 1]]; pf = -pf
        if A[k + 1, k] == 0: return 0.0
        pf *= A[k, k + 1]
        if k + 2 < n:
            tau = A[k, k + 2:] / A[k, k + 1]
            A[k + 2:, k + 2:] += np.outer(tau, A[k + 2:, k + 1]) - np.outer(A[k + 2:, k + 1], tau)
    return pf


def spin_qfi_profile(w, v, L, t, k=0, jmax=None):
    """Exact single-spin QFI at theta=0 for every read-out site j (0-based) after a rotation at site k.
       Returns (F_y[j], F_x[j]).  With sigma^x_j = prod_{l<j}(-i a_l b_l) a_j (and b_j for sigma^y_j),
           d<sigma^alpha_j>/dtheta = Re[-i s <0| sigma^alpha_j(t) g |0>],  s = (-1)^k, g = b_k (R_y) or a_k (R_x),
       and the vacuum expectation value of the Majorana string is a Pfaffian of two-point contractions.
       For j = 0 this reduces to F_y = G[a_1,a_k]^2 + G[b_1,a_k]^2; for j > 0 the string matters."""
    jmax = L if jmax is None else jmax
    G = Gt(w, v, t)
    S = np.zeros((2 * L, 2 * L)); idx = np.arange(L); S[idx, L + idx] = -1.0; S[L + idx, idx] = 1.0
    C0 = np.eye(2 * L) + 1j * S                      # <g_mu g_nu> in the vacuum (all spins down)
    Ctt = G @ C0 @ G.T                               # <g_mu(t) g_nu(t)>
    Ct0 = G @ C0                                     # <g_mu(t) g_nu(0)>
    order = np.empty(2 * L, int); order[0::2] = idx; order[1::2] = L + idx      # a_1 b_1 a_2 b_2 ...
    Fy = np.zeros(jmax); Fx = np.zeros(jmax); sgn = (-1.0) ** k
    for j in range(jmax):
        res = {}
        for enc, src in (('y', L + k), ('x', k)):
            dm2 = 0.0
            for last in (j, L + j):                  # sigma^x_j ends with a_j, sigma^y_j with b_j
                # <0| [string of 2j Majoranas at time t] g_last(t) g_src(0) |0> = Pf of the antisymmetrized matrix
                # of two-point contractions (Wick); the time-t operators use Ctt, the pair (t, 0) uses Ct0.
                ops = list(order[:2 * j]) + [last]
                n = len(ops) + 1
                M = np.zeros((n, n), complex)
                M[:n - 1, :n - 1] = Ctt[np.ix_(ops, ops)]
                M[:n - 1, n - 1] = Ct0[ops, src]
                M = np.triu(M, 1); M = M - M.T
                dm = np.real(-1j * sgn * (-1j) ** j * pfaffian(M))
                dm2 += dm * dm
            res[enc] = dm2
        Fy[j] = res['y']; Fx[j] = res['x']
    return Fy, Fx


def stationary_exact_ising(h, J=1.0):
    """gamma = 1, semi-infinite chain, T -> infinity: the left zero mode obeys phi_{j+1} = (h_j/J) phi_j, hence
       I_st(R_y, k=j=1) = [ sum_{j>=1} prod_{l<j} (h_l/J)^2 ]^(-2)   (zero if the sum diverges)."""
    lr = 2.0 * np.log(np.abs(np.asarray(h, float)) / J + 1e-300)
    c = np.concatenate(([0.0], np.cumsum(lr)[:-1]))
    m = c.max()
    return float(np.exp(-2.0 * (m + np.log(np.exp(c - m).sum()))))
