# FigR12: Matrix-product-state simulations: bond dimension and norm lost by truncation

FigR12 of the reply to the referees.

## Panels

(a) Bond dimension χ of |φ(t)⟩ (column chi_phi) versus tJ for L = 48 with δ = 0 (blue, to tJ = 15), 0.2J (orange, to tJ = 13) and 0.4J (green, to tJ = 3), and L = 100 with δ = 0.2J (black dashed, to tJ = 8). Dotted line: the maximum 128. At δ = 0.2J the curves for L = 48 and 100 coincide: chi_phi and chi_psi are equal at all 81 times tJ = 0, 0.1, …, 8; χ = 10 at tJ = 1, 99 at tJ = 6, and 128 from tJ = 6.8.

(b) Norm lost by truncation, 1 − ‖φ‖ (1 − norm_phi; norm_phi is the accumulated norm ‖φ‖ of the tangent state), of the same runs, log scale; points at tJ = 0 and with 1 − norm_phi ≤ 10⁻¹¹ are not drawn. At δ = 0.2J: 2.05 × 10⁻⁶ (L = 48) and 4.90 × 10⁻⁶ (L = 100) at tJ = 7; 6.36 × 10⁻⁶ and 1.62 × 10⁻⁵ at tJ = 8; 1.56 × 10⁻⁴ for L = 48 at tJ = 10.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| mps_L48_g1.00_h0.20_d0.00_chi128_dt0.1.txt | columns tJ, F_y, F_x, chi_psi, chi_phi, F_global, norm_phi; tJ = 0…15, step 0.1 | ../../codes/mps_edge_qfi.jl |
| mps_L48_g1.00_h0.20_d0.20_chi128_dt0.1.txt | same, tJ = 0…13 | ../../codes/mps_edge_qfi.jl |
| mps_L48_g1.00_h0.20_d0.40_chi128_dt0.1.txt | same, tJ = 0…3 (run stopped) | ../../codes/mps_edge_qfi.jl |
| mps_L100_g1.00_h0.20_d0.20_chi128_dt0.1.txt | same, tJ = 0…8 | ../../codes/mps_edge_qfi.jl |

## Reproduce

    python3 plot_FigR12.py

Needs numpy and matplotlib and the files in this folder.
