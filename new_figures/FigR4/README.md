# FigR4: Edge trace, averaging time, encoding-site scan and dependence on h

FigR4 of the reply to the referees.

## Panels

(a) Edge QFI F_y(t) (R_y) and F_x(t) (R_x) without averaging, k = 1, γ = 1, h = 0.5J, L = 400, 0 ≤ tJ ≤ 30 (stored to 60 with step 0.05). Dotted line: (1 - h²/J²)² = 0.5625. Inset: the same for h = 1.8J (trivial), log scale, tJ ≤ 6, dashed.

(b) I_1^st (solid) and the time-averaged QFI (dashed) for R_y versus T, 73 values T = 10^0 … 10^9 (logspace(0, 9, 73)), L = 8, 12, 16, 20, γ = 1, h = 0.5J, k = 1. Dotted lines: 0.5625 and 0.5625/2.

(c) I_1^st for R_y versus the encoding site k = 1…12, read-out at j = 1, L = 200, TJ = 10³, for (γ, h/J) = (1, 0.3), (1, 0.6), (0.5, 0.6). Lines: (1 - h²)² h^(2(k-1)) for γ = 1. Values below 3 × 10⁻⁷ are not drawn.

(d) I_1^st for R_y versus h/J = 0, 0.02, …, 1.4, k = 1, L = 400, TJ = 2 × 10³, γ = 1, 0.5, 0.25 (every second point drawn). Lines: stationary QFI of the clean semi-infinite chain (`free_fermion.plateau_exact`). Inset: points with h < J against 1 - h/J on log scales (every third point), with [8γ/(1+γ)²]² (1 - h/J)² (dashed) for 0.02 ≤ 1 - h/J ≤ 0.3.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| realtime_L400_g1.00_h0.50.txt | columns tJ, F_y, F_x; 1201 times 0…60; γ = 1, h = 0.5J, L = 400 | ../../codes/run_free.py (part A) |
| realtime_L400_g1.00_h1.80.txt | the same for h = 1.8J | ../../codes/run_free.py (part A) |
| window_scan_L{8,12,16,20}.txt | columns T, stat_y, stat_x, avg_y, avg_x; header gives the edge splitting ε₀ | ../../codes/run_free.py (part A) |
| encoding_scan.txt | columns γ, h, k, stat_y, stat_x, avg_y, avg_x, exact(γ = 1) | ../../codes/run_free.py (part B) |
| plateau_vs_h.txt | columns γ, h, stat_y, stat_x, avg_y, avg_x, exact_semi_infinite | ../../codes/run_free.py (part C) |

## Data

Generated from the repository root with:

    python3 codes/run_free.py

## Reproduce

    python3 plot_FigR4.py

Needs numpy and matplotlib and the files in this folder.
