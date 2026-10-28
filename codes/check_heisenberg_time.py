"""Time scales of the finite interacting chain that limit the lifetime fits: Heisenberg time and round-trip time.

Model as in all interacting runs (gamma = 1, Pauli matrices, open boundaries, J = 1):
    H = -(J/2) sum_j [(1+gamma) X_j X_{j+1} + (1-gamma) Y_j Y_{j+1}] - h sum_j Z_j + delta sum_j Z_j Z_{j+1}.

Heisenberg time. A finite chain has a discrete spectrum. For times shorter than t_H = 2 pi rho(E0), with rho the
density of many-body states at the energy E0 of the initial state, the levels are not resolved and the chain acts on
the edge mode like a continuum; the decay of the edge QFI is then the physical decay. Beyond t_H the evolution resolves
individual levels and the signal saturates at a finite-size value. Lifetimes must therefore be fitted at t << t_H.

Estimate used here (no diagonalization needed):
1. Even parity sector of H (it contains the initial state, all spins down), dimension M = 2^(L-1).
2. Exact first two moments of the spectrum in the sector: Ebar = Tr H / M, sigma^2 = Tr H^2 / M - Ebar^2
   (Tr H^2 is the sum of the squares of all matrix elements of the sparse H).
3. Gaussian density of states with these moments, rho(E) = M / (sqrt(2 pi) sigma) exp(-(E - Ebar)^2 / (2 sigma^2)),
   evaluated at E0 = <0|H|0>. The energy width of the initial state, sqrt(<0|H^2|0> - E0^2), is printed as well.
4. t_H = 2 pi rho(E0) (hbar = 1, J = 1).
This is an upper estimate: the reflection symmetry j -> L+1-j of the uniform chain splits the sector into two blocks of
about M/2 states, and only part of the states carries weight of the initial state. The true t_H may be about two times
shorter. The Gaussian form is accurate near the centre of the spectrum; E0 lies 1.5 to 2.6 sigma above it.

Round-trip time. The fastest single-particle excitations of the clean chain move with v = 2 min(h, J) sites per 1/J;
an excitation emitted at the edge returns after 2L/v.

Output: one line per (h, delta, L): h, delta, L, M, sigma, (E0 - Ebar)/sigma, energy width of |0>, rho(E0), t_H,
round-trip time 2L/v.
Run (from the folder of the revision):
    python3 codes/check_heisenberg_time.py > data/interacting/heisenberg_time.txt
Cost: builds the sparse sector Hamiltonian up to L = 22 (2^21 states); about 1 min, about 3 GB at L = 22.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_interacting as ed

print("# gamma=1, even sector; columns: h delta L dim sigma_H (E0-Ebar)/sigma width_of_|0> rho(E0) t_H round_trip")
for h, d in ((0.2, 0.2), (0.2, 0.3), (0.4, 0.1), (0.4, 0.2)):
    for L in (16, 18, 20, 22):
        st, H = ed.build_H_sector(L, 1.0, 1.0, h, d, 0); M = H.shape[0]
        Eb = H.diagonal().sum() / M                                   # Tr H / M
        sig = np.sqrt(H.multiply(H).sum() / M - Eb ** 2)              # sqrt(Tr H^2 / M - Ebar^2)
        v = np.zeros(M); v[np.searchsorted(st, 0)] = 1.0              # |0> = all spins down (bit string 0)
        Hv = H @ v; E0 = v @ Hv; w = np.sqrt(Hv @ Hv - E0 ** 2)      # energy and energy width of |0>
        rho = M / (np.sqrt(2 * np.pi) * sig) * np.exp(-(E0 - Eb) ** 2 / (2 * sig ** 2))
        print(f"{h:.2f} {d:.2f} {L} {M} {sig:.3f} {(E0 - Eb) / sig:+.3f} {w:.3f} {rho:.4g} {2 * np.pi * rho:.4g} "
              f"{2 * L / (2 * min(h, 1.0)):.1f}", flush=True)
