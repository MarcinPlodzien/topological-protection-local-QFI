# FigR3: Read-out scan of the submitted Fig. 2: zero-mode term and bulk wave packet

FigR3 of the reply to the referees.

## Panels

(a) Mean over tJ ∈ [150, 200] versus read-out site j = 1…60 (stored to 120) for L = 100 (solid), 200 (dashed), 400 (dotted); the zero-mode term (U_j0 U_10)², with 0 the mode of smallest energy, for L = 400 (dash-dotted); |mean − zero-mode term| for L = 400 (open circles).

(b) Excess of the mean over the zero-mode term, L = 400, windows tJ ∈ [50, 100], [150, 200], [450, 500]; values below 10⁻¹² are not drawn. Dotted vertical lines: j = 0.2 tJ at the end of each window (distance covered at the maximal group velocity 2 min(h, J) = 0.2J).

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| R1_data.txt | 120 rows (j); columns j, avg_L100, zm_L100, avg_L200, zm_L200, avg_L400, zm_L400 (avg: mean over [150, 200]; zm: zero-mode term; NaN for j > L) | ../../codes/compute_R1.py |
| R1_windows.txt | 120 rows (j); columns j, excess_50_100, excess_150_200, excess_450_500 (L = 400) | ../../codes/compute_R1.py |

## Data

Generated from the repository root with:

    python3 codes/compute_R1.py

## Reproduce

    python3 plot_FigR3.py

Needs numpy and matplotlib and the files in this folder.
