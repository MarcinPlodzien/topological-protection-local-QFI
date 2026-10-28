# FigR8: Topological island at the end of a trivial chain

FigR8 of the reply to the referees.

## Panels

(a) Edge-mode splitting ε₀ (smallest singular value of the chiral block B) versus island length n = 2…30. Line: ε₀(n = 12) × 0.5^(n-12), i.e. ∝ (h/J)^n with h = 0.5J, anchored at the row n = 12.

(b) Stationary QFI (QFI of the state averaged over [0, T]) versus T = 10^0.5 … 10^10 (77 values), R_y solid, R_x dashed, for n = 4, 8, 16, for the linear ramp (h_j rising linearly from 0.3J at j = 1 to 1.8J at j = 30, 1.8J beyond) and for the uniform topological chain (n = L, h = 0.5J).

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| island_splitting.txt | columns n, ε₀, stat_y(T = 30), stat_x(T = 30) | ../../codes/run_free.py (part D) |
| cases_island_n4.txt, cases_island_n8.txt, cases_island_n16.txt, cases_smooth_ramp.txt, cases_topological.txt | columns T, stat_y, stat_x, avg_y, avg_x; header: three lowest single-particle energies | ../../codes/run_free.py (part D) |

## Data

Generated from the repository root with:

    python3 codes/run_free.py

## Reproduce

    python3 plot_FigR8.py

Needs numpy and matplotlib and the files in this folder.
