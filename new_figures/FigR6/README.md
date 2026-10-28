# FigR6: Random fields and trivial edge states

FigR6 of the reply to the referees.

## Panels

(a) Random fields h_j = h + uniform[-Δ, Δ], h = 0.5J, L = 200, TJ = 10³, Δ/J = 0, 0.25, …, 4 (17 values), 1000 realizations per Δ (one at Δ = 0). Filled symbols: mean I_1^st for R_y (circles) and R_x (squares). Open symbols with dashed lines: mean time-averaged QFI for R_y and R_x. Band: interquartile range of I_1^st (R_y). Black line: exact mean of I_1^st for the semi-infinite chain at T → ∞, I_1^st = [Σ_j Π_{l<j} (h_l/J)²]⁻², 2 × 10⁴ sites, 40000 realizations per Δ, Δ = 0, 0.125, …, 4. Dashed vertical line and grey region: Δ > Δ_c = 2.6708J, the root of ⟨ln|h_j|⟩ = 0 (computed in the plot script by Delta_c() on a grid of 4 × 10⁵ points).

(b) I_1^st versus averaging time, R_y solid and R_x dashed, 77 values T = 10^0.5 … 10^10, L = 200: topological chain (h_j = 0.5J), trivial chain h_j = 1.8J with h_1 = 0.3J and with h_1 = 0, and islands whose first n = 8 or 16 sites have h = 0.5J and the rest 1.8J.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| disorder_L200_g1.00_h0.50_T1e+03.txt | 17 rows (Δ). Columns: 0 Δ; 1 mean Lyapunov exponent; 2–5 stat_y mean, median, q25, q75; 6–9 stat_x; 10–13 avg_y; 14–17 avg_x; 18 fraction with stat_y > 10 stat_x | ../../codes/run_disorder.py |
| semi_infinite_g1.00_h0.50.txt | 33 rows (Δ). Columns: Δ, mean, median, q25, q75 of the exact I_1^st, fraction with I_1^st > 10⁻⁶ | ../../codes/run_semi_infinite_refine.py |
| cases_topological.txt, cases_trivial_impurity.txt, cases_trivial_finetuned.txt, cases_island_n8.txt, cases_island_n16.txt | columns T, stat_y, stat_x, avg_y, avg_x; header: three lowest single-particle energies | ../../codes/run_free.py (part D) |

## Data

Generated from the repository root with:

    python3 codes/run_disorder.py 4
    python3 codes/run_semi_infinite_refine.py
    python3 codes/run_free.py

## Reproduce

    python3 plot_FigR6.py

Needs numpy and matplotlib and the files in this folder.
