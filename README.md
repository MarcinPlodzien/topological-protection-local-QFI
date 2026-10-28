# Topological protection of local quantum sensitivity: codes and data

M. Płodzień and J. Chwedeńczuk, "Topological protection of local quantum sensitivity", arXiv:2605.00770.

Codes, data, tests and plot scripts for the figures of the manuscript as submitted and for the new figures of the reply to the referees. Every figure and table has its own folder with the data it plots, a standalone plot script and the PDF it writes.

## Organization

    submitted_v1/   Figure1 ... Figure6 of the submitted manuscript: plotted data, plot script (make_*.py),
                    the PDF it writes, and the PDF as submitted (Figure<N>_as_submitted.pdf)
    new_figures/    new figures and tables of the reply to the referees, FigR1 ... FigR13, TableR1, TableR2,
                    each with data, plot_<name>.py, README.md and PDF; Fit_tau_delta/ (fit of the plateau lifetime)
    codes/          simulation and analysis codes (Python; one Julia script for matrix product states)
    data/           all data tables, plain text with a header line
      free/           free fermions (delta = 0): clean and engineered chains, single-spin maps
      disorder/       random fields: finite chains and exact semi-infinite results
      gap/            gap statistics of the disordered chain
      interacting/    exact propagation to tJ = 1000 (L = 12-22), full-diagonalization scans, autocorrelator,
                      lifetimes, level statistics, Heisenberg-time estimates;
                      trivial_edge/: chain with a trivial zero-energy edge state (FigR13)
      interacting_long/  exact propagation to tJ = 10^4 (L = 16, 18, 20)
      mps/            matrix-product-state results (L = 10, 48, 100)
      reply/          data of FigR3 and FigR5
    tests/          regression tests of the free-fermion and the interacting code

All paths in the scripts are relative to the script, so the folder can be moved as a whole.

## Reproduce the figures and tables

    cd new_figures/FigR1 && python3 plot_FigR1.py          # likewise FigR2 ... FigR13
    cd new_figures/TableR1 && python3 plot_TableR1.py      # prints Table R1
    cd new_figures/TableR2 && python3 plot_TableR2.py      # lifetimes from the traces, prints Table R2
    cd new_figures/Fit_tau_delta && python3 fit_tau_delta.py
    cd submitted_v1/Figure1 && python3 make_fig1_hero.py   # likewise Figure2 ... Figure6

Each plot script needs numpy and matplotlib and the files in its own folder, and runs in a few seconds. With matplotlib 3.10.8 the new figures FigR1 ... FigR13 and the submitted Figs. 1, 3, 4 and 6 are pixel-identical to the stored PDFs. For the submitted Figs. 2 and 5 the scripts plot the same data with different colors and legend placement.

## Regenerate the data

From the repository root. Scripts that write traces or tables skip existing files; delete a file to recompute it. Runtimes are given in the docstring of each script.

    # free fermions
    python3 codes/run_free.py
    python3 codes/run_phase_diagram.py
    python3 codes/run_spin_maps.py 4
    python3 codes/run_disorder.py 4
    python3 codes/run_semi_infinite_refine.py
    python3 codes/run_gap_disorder.py
    python3 codes/check_edge_cases.py > data/free/edge_cases.txt
    python3 codes/check_rx_finite_T.py > data/free/rx_finite_T.txt
    python3 codes/check_rx_decay.py > data/free/rx_decay_L1500.txt
    python3 codes/check_truncated_sum.py > data/disorder/truncated_sum_L200_h1.50.txt
    python3 codes/compute_R1.py                     # data of FigR3
    python3 codes/compute_R2.py                     # data of FigR5

    # interacting chain, exact propagation
    python3 codes/run_interacting.py 12             # L = 12-18 (and the gamma = 0.5, delta < 0 sets)
    python3 codes/run_interacting_long.py 20 2 6    # L = 20 to tJ = 1000
    python3 codes/run_interacting_long.py 22 2 6 8  # L = 22 to tJ = 1000
    python3 codes/run_long_traces.py 2 6            # L = 16, 18, 20 to tJ = 10^4
    python3 codes/run_xyz_line.py                   # integrable line h = 0, gamma = 0.5
    python3 codes/run_interacting_scan.py           # full-diagonalization scans in delta, L = 10, 12
    python3 codes/run_autocorrelator.py 12
    cd new_figures/FigR13 && python3 run_trivial_interacting.py 8 && python3 stationary_vs_T.py && cd ../..
    python3 codes/check_level_statistics.py 14 > data/interacting/level_statistics_L14.txt
    python3 codes/check_level_statistics.py 16 > data/interacting/level_statistics_L16_part1.txt
    python3 codes/check_heisenberg_time.py > data/interacting/heisenberg_time.txt
    python3 codes/check_submitted_fig6_L12.py > data/interacting/submitted_fig6_rescaled_L12.txt
    python3 codes/check_early_rates.py > data/interacting/early_rates.txt
    python3 codes/check_efold_times.py > data/interacting_long/efold_times.txt

    # matrix product states (Julia, ITensors/ITensorMPS)
    julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.2 15.0 0.1 128 1e-9
    julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.0 15.0 0.1 128 1e-9
    julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.4 15.0 0.1 128 1e-9
    julia codes/mps_edge_qfi.jl 100 1.0 0.2 0.2 8.0 0.1 128 1e-9
    julia codes/mps_edge_qfi.jl 10 1.0 0.2 0.2 6.0 0.05 64 1e-12   # benchmark against full diagonalization
    python3 codes/check_mps_chi_growth.py > data/mps/chi_growth.txt

    # all figures at once and the lifetime table
    python3 codes/make_figures.py                   # writes figures/FigR*.pdf and data/interacting/lifetimes.txt

## Model and conventions

H = -(J/2) Σ_j [(1+γ) X_j X_{j+1} + (1-γ) Y_j Y_{j+1}] - Σ_j h_j Z_j + δ Σ_j Z_j Z_{j+1}, with Pauli matrices X, Y, Z, open boundaries and J = 1 in all runs. The initial state has all spins down. The parameter θ is encoded by the rotation exp(-iθσ^a_k/2) of spin k (a = y for R_y, a = x for R_x; k = 1 unless stated), and the QFI of one spin j (j = 1 unless stated) is evaluated at θ = 0. F_y and F_x denote the edge QFI for R_y and R_x.

Free fermions (δ = 0): Jordan-Wigner Majoranas a_j = (Π_{l<j} Z_l) X_j and b_j = (Π_{l<j} Z_l) Y_j in the order γ = (a_1, …, a_L, b_1, …, b_L), so that H = (i/4) γ^T A γ with A real antisymmetric and γ(t) = G(t) γ(0), G(t) = exp(At). In this order A = [[0, B], [-B^T, 0]]; the SVD B = U diag(s) V^T gives the single-particle energies s_n and all blocks of G(t). For the edge spin F_y(t) = G[a_1,a_k]² + G[b_1,a_k]² and F_x(t) = G[a_1,b_k]² + G[b_1,b_k]².

Late-time quantities for a window [0, T]: "stat" (stationary QFI, I^st) is the QFI of the state averaged over [0, T]; "avg" is the QFI averaged over [0, T]. Both are evaluated in closed form with the kernel (1/T)∫_0^T e^{-ixt} dt; no time grid is involved.

Codes use 0-based site indices (k = 0 is spin 1); the Julia script uses 1-based indices. In the interacting code bit i of a basis index is the spin at site i, bit 1 = up.

## Software

Versions installed on the machine where the data were produced and checked (Intel Core i7-11800H, 16 logical CPUs, 62 GB RAM, Linux):

    Python 3.13.5 (Anaconda), numpy 2.1.3, scipy 1.15.3, matplotlib 3.10.8, numba 0.61.0
    Julia 1.12.5 with ITensors 0.6.23, ITensorMPS 0.2.6 (NDTensors 0.3.74)   -- only for codes/mps_edge_qfi.jl
    pdftoppm 22.02.0 (poppler)                                            -- only for the figure comparison below

The Python scripts need numpy and matplotlib; scipy is used by ed_interacting.py (sparse matrices, Bessel functions, eigsh), check_rx_decay.py and the tests; numba only by ed_interacting.cheb_traces_parity_fast (run_interacting_long.py, run_long_traces.py, test T9). The figure PDFs depend on the matplotlib version and on the serif font it finds; with matplotlib 3.10.8 they are pixel-identical to the PDFs in the folders.


## Tests

    cd tests
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 test_free.py          # 3.0 s, 0.15 GB
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 test_interacting.py   # 62 s, 0.78 GB (T9 needs numba)

test_free.py (codes/free_fermion.py, codes/bdg.py):

    F1  window kernel against trapezoidal integration of exp(-ixt) over [0, T]
    F2  chiral SVD fast path (late_time) against the 2L x 2L implementation (late_time_general), k = 1 and 4;
        stationary and time-averaged QFI against explicit time integration (L = 30, random fields, gamma = 0.7)
    F3  clean gamma = 1, L = 200: stationary R_y QFI = (1-h^2)^2, stationary R_x QFI = 0, h = 0.2, 0.5, 0.8
    F4  encoding-site scan: stationary QFI = (1-h^2)^2 h^(2(k-1)) for k = 2, 3, 5
    F5  zero-mode Lyapunov exponent: ln(h/J) for gamma = 1; exact roots of (1+g)x^2 - 2hx + (1-g) = 0 for gamma = 0.5, 0.3
    F6  closed form [sum_j prod_{l<j}(h_l/J)^2]^-2 against the numerical stationary QFI of a disordered chain; clean limit 0.5625
    F7  exact single-spin QFI at every site (Pfaffian formula) against exact diagonalization of the spin chain, L = 8

test_interacting.py (codes/ed_interacting.py):

    T1  bit-coded Hamiltonian against an independent construction from Kronecker products, L = 6
    T2  delta = 0: edge QFI of the interacting code against the free-fermion propagator of bdg.py, L = 8
    T3  delta = 0, gamma = 1: window average over 25 <= tJ <= 50 against (1-h^2/J^2)^2, L = 10
    T4  gamma = 1, h = 0: plateau against [1-(delta/J)^2]^2, L = 12
    T5  Chebyshev propagator against full diagonalization, norm conservation, L = 10
    T6  the submitted Trotter angles generate exchange J/2 and interaction delta/4: with these couplings the
        submitted Fig. 6 values are reproduced within 0.07 at L = 10, with the intended couplings they are missed by more than 0.1
    T7  parity-sector Chebyshev propagator against the full-space one, k = 3, j = 4, gamma = 0.7, L = 10
    T8  exact mapping (h, delta) -> (-h, -delta) of the edge QFI, L = 10
    T9  multithreaded propagator cheb_traces_parity_fast (numba, 2 threads) against cheb_traces_parity,
        gamma = 0.7, h = 0.3J, delta = 0.25J, L = 10, tJ <= 20

Result: every check prints PASS, and both files end with "ALL TESTS PASSED" (exit code 0). Largest deviations: F7 2.6e-15, T7 1.9e-16, T8 3.7e-14, T9 2.3e-15 (norm 1.0e-15); T6 prints 0.7030, 0.6402, 0.9008 for the rescaled couplings against the submitted 0.7070, 0.7040, 0.9058 (L = 10 against the submitted L = 12).

