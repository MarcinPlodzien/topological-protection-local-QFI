# FigR10: Interacting chain: size dependence, matrix product states, short-chain scans, sign of delta

FigR10 of the reply to the referees.

## Panels

(a) Running average of F_y(t), γ = 1, for L = 12, 14, 16, 18 (grey to black), 20 (magenta) and 22 (indigo), (h, δ) = (0.2, 0.3)J solid, (0.4, 0.2)J dashed, (0.2, 0.15)J dotted; 5 ≤ tJ ≤ 1000.

(b) MPS (two-site TDVP), L = 48, γ = 1, h = 0.2J, time step 0.1/J, bond dimension χ ≤ 128, cutoff 10⁻⁹. Lines: F_y (upper) and F_x (lower). δ = 0.2J (magenta) is drawn up to tJ = 10 and compared with exact propagation for L = 18 (circles, every 0.5/J); δ = 0 (indigo) to tJ = 15 and compared with the free-fermion result for L = 48 computed in the plot script from the SVD of B (squares, every 0.5/J).

(c) Window averages of F_y versus δ = -0.6 … 0.6 (step 0.025), γ = 1, h = 0.2J, full diagonalization, L = 10 dashed and L = 12 solid; windows tJ ∈ [10, 30], [30, 100], [100, 300], [300, 1000], each averaged over 120 random times drawn with np.random.default_rng(1). Dotted grey: [1 - (δ/J)²]², the exact plateau at h = 0.

(d) Running averages of F_y (solid) and F_x (dotted), L = 16: γ = 1, h = 0.2J, δ = ±0.2J; γ = 0.5, h = 0.2J, δ = 0.1J and 0.2J.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| trace_L{12,14,16,18}_g1.00_h0.20_d0.30.txt | columns tJ, F_y, F_x; step 0.5/J; to tJ = 1000 | ../../codes/run_interacting.py |
| trace_L{12,14,16,18}_g1.00_h0.40_d0.20.txt | same | ../../codes/run_interacting.py |
| trace_L{12,14,16,18}_g1.00_h0.20_d0.15.txt | same | ../../codes/run_interacting.py |
| trace_L{20,22}_g1.00_h0.20_d0.30.txt, trace_L{20,22}_g1.00_h0.40_d0.20.txt, trace_L{20,22}_g1.00_h0.20_d0.15.txt | same, to tJ = 1000 | ../../codes/run_interacting_long.py |
| trace_L18_g1.00_h0.20_d0.20.txt | same (panel b) | ../../codes/run_interacting.py |
| trace_L16_g1.00_h0.20_d0.20.txt, trace_L16_g1.00_h0.20_d-0.20.txt, trace_L16_g0.50_h0.20_d0.10.txt, trace_L16_g0.50_h0.20_d0.20.txt | same (panel d) | ../../codes/run_interacting.py |
| mps_L48_g1.00_h0.20_d0.20_chi128_dt0.1.txt | columns tJ, F_y, F_x, χ(ψ), χ(φ), F_global, norm_φ; tJ = 0…13 step 0.1 | ../../codes/mps_edge_qfi.jl |
| mps_L48_g1.00_h0.20_d0.00_chi128_dt0.1.txt | same, δ = 0, tJ = 0…15 | ../../codes/mps_edge_qfi.jl |
| scan_L10_g1.00_h0.20.txt, scan_L12_g1.00_h0.20.txt | 49 rows (δ); columns δ, mean F_y in the six windows [10,30], [30,100], [100,300], [300,10³], [10³,3·10³], [3·10³,10⁴], then mean F_x in the same windows | ../../codes/run_interacting_scan.py |

## Data

Generated from the repository root with:

    python3 codes/run_interacting.py 12
    python3 codes/run_interacting_long.py 20 2 6
    python3 codes/run_interacting_long.py 22 2 6 8
    python3 codes/run_interacting_scan.py

## Reproduce

    python3 plot_FigR10.py

Needs numpy and matplotlib and the files in this folder.
