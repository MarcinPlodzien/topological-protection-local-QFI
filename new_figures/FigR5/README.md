# FigR5: Edge QFI near the transition in the setting of the submitted Fig. 4

FigR5 of the reply to the referees.

## Panels

(a) Time-averaged edge QFI, the mean of F_y over the 100 times linspace(200, 300, 100) used by the submitted generator, versus 1 − h/J = 10⁻³ … 0.5 (50 log-spaced values), for L = 50, 100, 200, 400. Dotted: (1 − h²/J²)². Dashed vertical line: 1 − h/J = 0.02, the lower end of the submitted figure.

(b) F_y(t) at h = 0.98J for the four lengths, tJ = 0 … 450 (1801 times, step 0.25). Shaded: the window [200, 300]. Dotted: (1 − 0.98²)² = 1.58 × 10⁻³.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| R2_data.txt | 50 rows; columns one_minus_h, avg_L50, stat_L50, avg_L100, stat_L100, avg_L200, stat_L200, avg_L400, stat_L400, exact_semi_infinite (avg: mean over the 100 times; stat: QFI of the state averaged over [200, 300], not plotted) | ../../codes/compute_R2.py |
| R2_trace.txt | 1801 rows; columns tJ, F_L50, F_L100, F_L200, F_L400 at h = 0.98J | ../../codes/compute_R2.py |

## Data

Generated from the repository root with:

    python3 codes/compute_R2.py

## Reproduce

    python3 plot_FigR5.py

Needs numpy and matplotlib and the files in this folder.
