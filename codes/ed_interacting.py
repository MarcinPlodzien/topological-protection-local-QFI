#!/usr/bin/env python3
"""
ed_interacting.py -- single-site QFI in the interacting chain, numerically exact.

    H = -(J/2) sum_j [(1+gamma) sx_j sx_{j+1} + (1-gamma) sy_j sy_{j+1}]
        - h sum_j sz_j + delta sum_j sz_j sz_{j+1}          (Pauli matrices, OBC)

All coefficients multiply Pauli strings.  The state |dd...d> is rotated at site k
by R_a(theta) = exp(-i theta s^a_k / 2), a = y or x, evolved with H, and the QFI
for theta (at theta = 0) of the reduced state of site j is returned.

Four propagators:
  * dense_traces  : full diagonalization, any time grid (L <= 14)
  * cheb_traces   : Chebyshev expansion of exp(-iHt) with sparse H (L <= 22), full Hilbert space
  * cheb_traces_parity : the same in the two parity sectors (H conserves the parity of the number of up
                    spins; |psi0> is even, the flipped state is odd), about three times faster
  * cheb_traces_parity_fast : cheb_traces_parity with a numba kernel parallel over rows (needs numba);
                    used for L = 20, 22 to tJ = 1000 and for the traces to tJ = 10^4
All are checked against each other and against closed forms in tests/.

Basis: bit i of the basis index is the spin at site i (0-based); bit = 1 is up.

Method.  d|psi_theta>/dtheta at theta = 0 is proportional to |phi> = (spin k flipped)|psi0>, so the QFI needs
only the two propagated states U(t)|psi0> and U(t)|phi>: the reduced state of site j gives the Bloch vector m
and its derivative dm (from <psi|sigma_j|phi>), and the qubit QFI is |dm|^2 + (m.dm)^2/(1-|m|^2) (qfi_site).
No finite differences in theta are taken.

Units and conventions: J = 1 in all runs; sites 0-based (k = j = 0 is the edge spin of the paper); the initial
state is all spins down (basis index 0).  Used by run_interacting.py, run_interacting_scan.py, run_xyz_line.py,
run_autocorrelator.py, run_interacting_long.py, run_long_traces.py, check_submitted_fig6_L12.py,
check_level_statistics.py, check_heisenberg_time.py and tests/test_interacting.py.  No input or output files.
Memory: 2^L complex vectors (16 MB at L = 20) and a sparse H with L entries per row; dense_traces needs the
full 2^L x 2^L eigenvector matrix (128 MB at L = 12, 2 GB at L = 14).
Test T9: cheb_traces_parity_fast equals cheb_traces_parity to 2e-15 at L = 10.
"""
import os
import numpy as np
import scipy.sparse as sp
from scipy.special import jv
from scipy.sparse.linalg import eigsh


# ----------------------------------------------------------------------------
#  Hamiltonian
# ----------------------------------------------------------------------------
def build_H(L, J, gamma, h, delta):
    """Return (diag, offdiag) with H = diag(diag) + offdiag; offdiag is real CSR."""
    N = 1 << L
    idx = np.arange(N, dtype=np.int64)
    z = [1.0 - 2.0 * (1 - ((idx >> i) & 1)) for i in range(L)]      # +1 up, -1 down
    h = np.broadcast_to(np.asarray(h, float), (L,))
    diag = np.zeros(N)
    for i in range(L):
        diag -= h[i] * z[i]
    cols = np.empty((N, L - 1), dtype=np.int32)
    vals = np.empty((N, L - 1))
    for i in range(L - 1):
        diag += delta * z[i] * z[i + 1]
        differ = (z[i] * z[i + 1]) < 0
        cols[:, i] = idx ^ ((1 << i) | (1 << (i + 1)))
        # XX flips both spins with +1; YY flips both with +1 (antiparallel) or -1 (parallel)
        vals[:, i] = -0.5 * J * ((1 + gamma) + (1 - gamma) * np.where(differ, 1.0, -1.0))
    off = sp.csr_matrix((vals.ravel(), cols.ravel(), np.arange(0, N * (L - 1) + 1, L - 1)),
                        shape=(N, N))
    return diag, off


def initial_states(L, k):
    """|psi0> = all down; |phi> = spin k flipped.  d|psi>/dtheta = -phi/2 (R_y), -i phi/2 (R_x)."""
    N = 1 << L
    psi0 = np.zeros(N, complex); psi0[0] = 1.0
    phi = np.zeros(N, complex); phi[1 << k] = 1.0
    return psi0, phi


# ----------------------------------------------------------------------------
#  Single-site QFI at theta = 0
# ----------------------------------------------------------------------------
def _bloch(a, b, L, j):
    """m = <a|sigma_j|a>,  c = <a|sigma_j|b> (complex 3-vector)."""
    A = np.moveaxis(a.reshape((2,) * L), L - 1 - j, 0).reshape(2, -1)   # A[s, rest], s = bit j
    B = np.moveaxis(b.reshape((2,) * L), L - 1 - j, 0).reshape(2, -1)
    def vec(r):                       # r[s, s'] = <s| ... |s'> matrix elements in bit basis (1 = up)
        up, dn = 1, 0
        x = r[up, dn] + r[dn, up]
        y = -1j * r[up, dn] + 1j * r[dn, up]
        zc = r[up, up] - r[dn, dn]
        return np.array([x, y, zc])
    raa = A.conj() @ A.T              # raa[s, s'] = sum_rest conj(A[s]) A[s']
    rab = A.conj() @ B.T
    return vec(raa).real, vec(rab)


def qfi_site(psi, phi, L, j):
    """QFI of site j for R_y and R_x encoding; psi = U|psi0>, phi = U|flipped>."""
    m, c = _bloch(psi, phi, L, j)
    out = []
    # qubit QFI for a state with Bloch vector m and derivative dm: |dm|^2 + (m.dm)^2/(1 - |m|^2)
    for pref in (-0.5, -0.5j):        # d psi = pref * phi
        dm = 2.0 * np.real(pref * c)
        out.append(dm @ dm + (m @ dm) ** 2 / max(1.0 - m @ m, 1e-14))
    return out                        # [F_y, F_x]


# ----------------------------------------------------------------------------
#  Propagators
# ----------------------------------------------------------------------------
def dense_traces(L, J, gamma, h, delta, times, k=0, j=0):
    diag, off = build_H(L, J, gamma, h, delta)
    E, V = np.linalg.eigh(off.toarray() + np.diag(diag))
    psi0, phi = initial_states(L, k)
    c0 = V.T @ psi0; c1 = V.T @ phi
    out = np.empty((len(times), 2))
    for n, t in enumerate(times):
        ph = np.exp(-1j * E * t)
        out[n] = qfi_site(V @ (ph * c0), V @ (ph * c1), L, j)
    return out


def cheb_traces(L, J, gamma, h, delta, tmax, dt=0.5, nsub=4, k=0, j=0, tol=1e-12, progress=None):
    """Times dt, 2dt, ..., tmax.  One Chebyshev expansion per block of nsub outputs."""
    diag, off = build_H(L, J, gamma, h, delta)
    Hop = sp.linalg.LinearOperator((1 << L, 1 << L), matvec=lambda v: diag * v + off @ v, dtype=float)
    emax = eigsh(Hop, k=1, which='LA', tol=1e-3, return_eigenvectors=False)[0]
    emin = eigsh(Hop, k=1, which='SA', tol=1e-3, return_eigenvectors=False)[0]
    # spectrum mapped into [-1, 1]: H = a Ht + b, with 2 % margin on the extremal eigenvalues
    a = 0.5 * (emax - emin) * 1.02; b = 0.5 * (emax + emin)
    def Ht(v):
        return (diag * v + off @ v - b * v) / a
    # Chebyshev expansion  exp(-i H tau) = exp(-i b tau) sum_k (2 - delta_k0) (-i)^k J_k(a tau) T_k(Ht);
    # the Bessel functions J_k(a tau) decay faster than exponentially for k > a tau, and the series is cut
    # where |J_kmax(a tau)| < tol (1e-12).  One set of T_k(Ht) v serves the nsub output times tau = dt ... nsub dt.
    taus = dt * np.arange(1, nsub + 1)
    kmax = int(a * taus[-1] + 12 * (a * taus[-1]) ** (1 / 3) + 25)
    while abs(jv(kmax, a * taus[-1])) > tol:
        kmax += 5
    ks = np.arange(kmax + 1)
    coef = np.array([(2.0 - (ks == 0)) * (-1j) ** ks * jv(ks, a * tau) * np.exp(-1j * b * tau)
                     for tau in taus])                                   # (nsub, kmax+1)
    def block(v):
        outs = [coef[s, 0] * v for s in range(nsub)]
        t0, t1 = v, Ht(v)
        for s in range(nsub):
            outs[s] += coef[s, 1] * t1
        for n in range(2, kmax + 1):                     # recursion T_n = 2 Ht T_{n-1} - T_{n-2}
            t0, t1 = t1, 2.0 * Ht(t1) - t0
            for s in range(nsub):
                if abs(coef[s, n]) > 1e-17:
                    outs[s] += coef[s, n] * t1
        return outs
    psi, phi = initial_states(L, k)
    nblocks = int(round(tmax / taus[-1]))
    times = []; res = []
    for nb in range(nblocks):
        P = block(psi); Q = block(phi)
        for s in range(nsub):
            times.append(nb * taus[-1] + taus[s]); res.append(qfi_site(P[s], Q[s], L, j))
        psi, phi = P[-1], Q[-1]
        if progress and (nb + 1) % progress == 0:
            print(f"   L={L} h={h} d={delta}: t={times[-1]:.0f} norm={np.linalg.norm(psi):.10f}", flush=True)
    return np.array(times), np.array(res), float(np.linalg.norm(psi))


def build_H_sector(L, J, gamma, h, delta, parity):
    """H restricted to the states with an even (parity=0) or odd (1) number of up spins.
       Returns (states, H) with H real CSR including the diagonal."""
    N = 1 << L
    idx = np.arange(N, dtype=np.int64)
    par = np.zeros(N, dtype=np.int8)
    for i in range(L):
        par ^= ((idx >> i) & 1).astype(np.int8)
    st = idx[par == parity]; M = len(st)
    pos = np.empty(N, dtype=np.int32); pos[st] = np.arange(M, dtype=np.int32)
    z = [2.0 * ((st >> i) & 1) - 1.0 for i in range(L)]
    h = np.broadcast_to(np.asarray(h, float), (L,))
    diag = np.zeros(M)
    for i in range(L):
        diag -= h[i] * z[i]
    cols = np.empty((M, L), dtype=np.int32); vals = np.empty((M, L))
    for i in range(L - 1):
        diag += delta * z[i] * z[i + 1]
        cols[:, i] = pos[st ^ ((1 << i) | (1 << (i + 1)))]
        vals[:, i] = -0.5 * J * ((1 + gamma) + (1 - gamma) * np.where(z[i] * z[i + 1] < 0, 1.0, -1.0))
    cols[:, L - 1] = np.arange(M); vals[:, L - 1] = diag
    H = sp.csr_matrix((vals.ravel(), cols.ravel(), np.arange(0, M * L + 1, L)), shape=(M, M))
    H.sort_indices()
    return st, H


def cheb_traces_parity(L, J, gamma, h, delta, tmax, dt=0.5, nsub=4, k=0, j=0, tol=1e-12, progress=None,
                       checkpoint=None, every=25):
    """Same output as cheb_traces.  |psi0> is propagated in the even sector, the flipped state in the odd one.
    checkpoint: path of an .npz file; the state and the trace are saved there after every `every` blocks
    (every*nsub*dt units of time) and the run continues from it if the file exists."""
    # H conserves the parity of the number of up spins (all terms flip spins in pairs or are diagonal):
    # |psi0> (all down) lies in the even sector, the flipped state |phi> in the odd one.  Each sector has
    # dimension 2^(L-1); the common spectral bounds of both sectors set the Chebyshev scaling.
    sec = [build_H_sector(L, J, gamma, h, delta, p) for p in (0, 1)]
    lo = min(eigsh(H, k=1, which='SA', tol=1e-3, return_eigenvectors=False)[0] for _, H in sec)
    hi = max(eigsh(H, k=1, which='LA', tol=1e-3, return_eigenvectors=False)[0] for _, H in sec)
    a = 0.5 * (hi - lo) * 1.02; b = 0.5 * (hi + lo)
    Hs = [((H - b * sp.identity(H.shape[0], format='csr')) * (2.0 / a)).tocsr() for _, H in sec]   # 2 (H - b)/a
    taus = dt * np.arange(1, nsub + 1)
    kmax = int(a * taus[-1] + 12 * (a * taus[-1]) ** (1 / 3) + 25)
    while abs(jv(kmax, a * taus[-1])) > tol:
        kmax += 5
    ks = np.arange(kmax + 1)
    coef = np.array([(2.0 - (ks == 0)) * (-1j) ** ks * jv(ks, a * tau) * np.exp(-1j * b * tau) for tau in taus])
    nmax = [int(np.nonzero(np.abs(c) > 1e-17)[0].max()) for c in coef]
    def block(H2, v):
        # H2 = 2 (H - b)/a, so T_1 v = (H2/2) v and T_n = H2 T_{n-1} - T_{n-2}, in place to save memory
        outs = [coef[s, 0] * v for s in range(nsub)]
        t0 = v; t1 = 0.5 * (H2 @ v)
        for s in range(nsub):
            outs[s] += coef[s, 1] * t1
        for n in range(2, kmax + 1):
            t2 = H2 @ t1; t2 -= t0
            t0, t1 = t1, t2
            for s in range(nsub):
                if n <= nmax[s]:
                    outs[s] += coef[s, n] * t1
        return outs
    N = 1 << L
    psi = np.zeros(len(sec[0][0]), complex); psi[np.searchsorted(sec[0][0], 0)] = 1.0
    phi = np.zeros(len(sec[1][0]), complex); phi[np.searchsorted(sec[1][0], 1 << k)] = 1.0
    full_a = np.zeros(N, complex); full_b = np.zeros(N, complex)
    nblocks = int(round(tmax / taus[-1])); times = []; res = []; nb0 = 0
    if checkpoint and os.path.exists(checkpoint):
        c = np.load(checkpoint)
        psi, phi, nb0 = c['psi'], c['phi'], int(c['nb'])
        times, res = list(c['times']), [r for r in c['res']]
    for nb in range(nb0, nblocks):
        P = block(Hs[0], psi); Q = block(Hs[1], phi)
        for s in range(nsub):
            # embed the sector vectors into the full 2^L space for the reduced state of site j
            full_a[sec[0][0]] = P[s]; full_b[sec[1][0]] = Q[s]
            times.append(nb * taus[-1] + taus[s]); res.append(qfi_site(full_a, full_b, L, j))
        psi, phi = P[-1], Q[-1]
        if progress and (nb + 1) % progress == 0:
            print(f"   L={L} h={h} d={delta}: t={times[-1]:.0f} norm={np.linalg.norm(psi):.10f}", flush=True)
        if checkpoint and (nb + 1) % every == 0 and nb + 1 < nblocks:
            np.savez(checkpoint + '.tmp.npz', psi=psi, phi=phi, nb=nb + 1, times=np.array(times), res=np.array(res))
            os.replace(checkpoint + '.tmp.npz', checkpoint)
            np.savetxt(checkpoint[:-4] + '.partial.txt', np.column_stack([times, res]), fmt="%.6f %.8e %.8e",
                       header=f"PARTIAL trace, L={L} gamma={gamma} h={h} delta={delta}; tJ  F_y  F_x")
    if checkpoint:
        for f in (checkpoint, checkpoint[:-4] + '.partial.txt'):
            if os.path.exists(f):
                os.remove(f)
    return np.array(times), np.array(res), float(np.linalg.norm(psi))


# ---------------------------------------------------------------------------------------------------------------
# Multithreaded variant of cheb_traces_parity.
# Same algorithm, same Chebyshev coefficients and truncation, same output. Differences are only in the
# implementation: the sparse product H2 @ t1, the recurrence t2 = H2 t1 - t0 and the accumulation of the nsub
# outputs are fused in one numba kernel that runs over the rows in parallel, works directly on complex vectors
# (scipy converts the real matrix to complex on every product) and reuses three buffers.
try:
    import numba as _nb

    @_nb.njit(parallel=True, cache=True)
    def _cheb_step(indptr, indices, data, t1, t0, out, outs, c):
        # out[i] = sum_k data[k] t1[indices[k]] - t0[i];  outs[s, i] += c[s] * out[i]
        M = t1.shape[0]; ns = outs.shape[0]
        for i in _nb.prange(M):
            acc = 0j
            for p in range(indptr[i], indptr[i + 1]):
                acc += data[p] * t1[indices[p]]
            acc -= t0[i]
            out[i] = acc
            for s in range(ns):
                outs[s, i] += c[s] * acc

    @_nb.njit(parallel=True, cache=True)
    def _half_matvec(indptr, indices, data, v, out):
        for i in _nb.prange(v.shape[0]):
            acc = 0j
            for p in range(indptr[i], indptr[i + 1]):
                acc += data[p] * v[indices[p]]
            out[i] = 0.5 * acc
except ImportError:                     # numba missing: cheb_traces_parity_fast is unavailable
    _nb = None


def cheb_traces_parity_fast(L, J, gamma, h, delta, tmax, dt=0.5, nsub=4, k=0, j=0, tol=1e-12, progress=None,
                            checkpoint=None, every=25, threads=None):
    """Same output as cheb_traces_parity (checked against it in tests/test_interacting.py, T9).
    threads: number of numba threads (default: numba's default). checkpoint/every as in cheb_traces_parity."""
    if _nb is None:
        raise ImportError("cheb_traces_parity_fast needs numba")
    if threads:
        _nb.set_num_threads(int(threads))
    sec = [build_H_sector(L, J, gamma, h, delta, p) for p in (0, 1)]
    lo = min(eigsh(H, k=1, which='SA', tol=1e-3, return_eigenvectors=False)[0] for _, H in sec)
    hi = max(eigsh(H, k=1, which='LA', tol=1e-3, return_eigenvectors=False)[0] for _, H in sec)
    a = 0.5 * (hi - lo) * 1.02; b = 0.5 * (hi + lo)
    Hs = [((H - b * sp.identity(H.shape[0], format='csr')) * (2.0 / a)).tocsr() for _, H in sec]   # 2 (H - b)/a
    for H2 in Hs:
        H2.sort_indices()
    taus = dt * np.arange(1, nsub + 1)
    kmax = int(a * taus[-1] + 12 * (a * taus[-1]) ** (1 / 3) + 25)
    while abs(jv(kmax, a * taus[-1])) > tol:
        kmax += 5
    ks = np.arange(kmax + 1)
    coef = np.array([(2.0 - (ks == 0)) * (-1j) ** ks * jv(ks, a * tau) * np.exp(-1j * b * tau) for tau in taus])
    nmax = [int(np.nonzero(np.abs(c) > 1e-17)[0].max()) for c in coef]
    # coefficients beyond nmax[s] are dropped exactly as in cheb_traces_parity
    C = np.array([[coef[s, n] if n <= nmax[s] else 0.0 for s in range(nsub)] for n in range(kmax + 1)], complex)

    arrs = [(H2.indptr.astype(np.int64), H2.indices.astype(np.int64), H2.data) for H2 in Hs]   # converted once

    def block(sector, v):
        ip, ix, dd = arrs[sector]
        M = v.shape[0]
        outs = np.empty((nsub, M), complex)
        for s in range(nsub):
            outs[s] = C[0, s] * v
        t0 = v.copy(); t1 = np.empty(M, complex); t2 = np.empty(M, complex)
        _half_matvec(ip, ix, dd, v, t1)                  # t1 = 0.5 H2 v
        for s in range(nsub):
            outs[s] += C[1, s] * t1
        for n in range(2, kmax + 1):
            _cheb_step(ip, ix, dd, t1, t0, t2, outs, C[n])
            t0, t1, t2 = t1, t2, t0
        return list(outs)

    N = 1 << L
    psi = np.zeros(len(sec[0][0]), complex); psi[np.searchsorted(sec[0][0], 0)] = 1.0
    phi = np.zeros(len(sec[1][0]), complex); phi[np.searchsorted(sec[1][0], 1 << k)] = 1.0
    full_a = np.zeros(N, complex); full_b = np.zeros(N, complex)
    nblocks = int(round(tmax / taus[-1])); times = []; res = []; nb0 = 0
    if checkpoint and os.path.exists(checkpoint):
        c = np.load(checkpoint)
        psi, phi, nb0 = c['psi'], c['phi'], int(c['nb'])
        times, res = list(c['times']), [r for r in c['res']]
    for nb in range(nb0, nblocks):
        P = block(0, psi); Q = block(1, phi)
        for s in range(nsub):
            full_a[sec[0][0]] = P[s]; full_b[sec[1][0]] = Q[s]
            times.append(nb * taus[-1] + taus[s]); res.append(qfi_site(full_a, full_b, L, j))
        psi, phi = P[-1], Q[-1]
        if progress and (nb + 1) % progress == 0:
            print(f"   L={L} h={h} d={delta}: t={times[-1]:.0f} norm={np.linalg.norm(psi):.10f}", flush=True)
        if checkpoint and (nb + 1) % every == 0 and nb + 1 < nblocks:
            np.savez(checkpoint + '.tmp.npz', psi=psi, phi=phi, nb=nb + 1, times=np.array(times), res=np.array(res))
            os.replace(checkpoint + '.tmp.npz', checkpoint)
            np.savetxt(checkpoint[:-4] + '.partial.txt', np.column_stack([times, res]), fmt="%.6f %.8e %.8e",
                       header=f"PARTIAL trace, L={L} gamma={gamma} h={h} delta={delta}; tJ  F_y  F_x")
    if checkpoint:
        for f in (checkpoint, checkpoint[:-4] + '.partial.txt'):
            if os.path.exists(f):
                os.remove(f)
    return np.array(times), np.array(res), float(np.linalg.norm(psi))
