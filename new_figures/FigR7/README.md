# FigR7: Disordered chain: mean, median, dependence on L and T, gamma = 0.5, trivial chain

FigR7 of the reply to the referees.

## Panels

(a) γ = 1, h = 0.5J. Mean (circles) and median (squares) of the stationary QFI I_1^st (R_y) for L = 200, TJ = 10³. Lines: mean and median of the exact I_1^st = [Σ_j Π_{l<j}(h_l/J)²]⁻² of the semi-infinite chain at T → ∞ (2 × 10⁴ sites, 40000 realizations per Δ, Δ = 0, 0.125, …, 4); band: its interquartile range. Dashed line: Δ_c = 2.6708J.

(b) Mean I_1^st (R_y), γ = 1, h = 0.5J, for L = 100, 200, 400 at TJ = 10³ and for L = 200 at TJ = 10² and 10⁴, with the exact line.

(c) γ = 0.5, h = 0.5J, L = 200, TJ = 10³: mean stationary QFI (filled) and mean time-averaged QFI (open) for R_y and R_x. Dashed line: zero of the disorder-averaged Lyapunov exponent of the zero mode (column 1), by linear interpolation between grid points.

(d) γ = 1, h = 1.5J, L = 200, TJ = 10³, same quantities. Shaded: 1.700 < Δ/J < 1.917, where ⟨ln|h_j|⟩ < 0 and the random chain is topological. Black line: exact mean for the semi-infinite chain (4000 realizations, Δ = 1, 1.05, …, 3).

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| disorder_L200_g1.00_h0.50_T1e+03.txt | 17 rows (Δ); columns 0 Δ; 1 mean Lyapunov exponent; 2–5 stat_y mean, median, q25, q75; 6–9 stat_x; 10–13 avg_y; 14–17 avg_x; 18 fraction with stat_y > 10 stat_x | ../../codes/run_disorder.py |
| disorder_L100_g1.00_h0.50_T1e+03.txt | same, L = 100 | ../../codes/run_disorder.py |
| disorder_L400_g1.00_h0.50_T1e+03.txt | same, L = 400, 300 realizations | ../../codes/run_disorder.py |
| disorder_L200_g1.00_h0.50_T1e+02.txt, disorder_L200_g1.00_h0.50_T1e+04.txt | same, TJ = 10², 10⁴ | ../../codes/run_disorder.py |
| disorder_L200_g0.50_h0.50_T1e+03.txt | same, γ = 0.5 | ../../codes/run_disorder.py |
| disorder_L200_g1.00_h1.50_T1e+03.txt | same, h = 1.5J | ../../codes/run_disorder.py |
| semi_infinite_g1.00_h0.50.txt | 33 rows (Δ); columns Δ, mean, median, q25, q75 of the exact I_1^st, fraction with I_1^st > 10⁻⁶; 40000 realizations | ../../codes/run_semi_infinite_refine.py |
| semi_infinite_g1.00_h1.50.txt | 41 rows (Δ = 1…3); same columns; 4000 realizations | ../../codes/run_disorder.py |

## Data

Generated from the repository root with:

    python3 codes/run_disorder.py 4
    python3 codes/run_semi_infinite_refine.py

## Reproduce

    python3 plot_FigR7.py

Needs numpy and matplotlib and the files in this folder.
