# FigR1: Phase diagram of the stationary edge QFI and exact single-spin QFI maps

FigR1 of the reply to the referees.

## Panels

(a) Stationary edge QFI I_1^st (QFI of the state of spin 1 averaged over the window [0, T]) for R_y and read-out j = 1, over the grid h/J = 0, 0.025, …, 2 (81 values) and γ = 0, 0.025, …, 1 (41 values). L = 200, TJ = 200. Colour scale linear from 0 to 1. The dashed line is h = J; the labels w = 1 and w = 0 mark the topological and trivial phase.

(b) Exact QFI I_j(t) of spin j (Jordan-Wigner string included, Pfaffian formula) after the rotation of spin 1, for j = 1…60 and 121 times tJ = 0, 0.5, …, 60. L = 100, γ = 0.3, h = 0.5J (topological). The x axis is cut at j = 40.

(c) The same for γ = 0.8, h = 1.8J (trivial).

Panels (b) and (c) share one square-root colour scale (PowerNorm(0.5)) from 0 to the larger maximum of the two maps.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| phase_diagram_h.txt | the 81 values of h/J (columns of the table below) | ../../codes/run_phase_diagram.py |
| phase_diagram_gamma.txt | the 41 values of γ (rows) | ../../codes/run_phase_diagram.py |
| phase_diagram_stat_y.txt | 41 × 81 table of I_1^st for R_y, rows γ, columns h | ../../codes/run_phase_diagram.py |
| spinmap_times.txt | the 121 times tJ | ../../codes/run_spin_maps.py |
| spinmap_topological_Ry.txt | 121 × 60 table, rows time, columns j = 1…60: single-spin QFI, γ = 0.3, h = 0.5J | ../../codes/run_spin_maps.py |
| spinmap_trivial_Ry.txt | the same for γ = 0.8, h = 1.8J | ../../codes/run_spin_maps.py |

## Data

Generated from the repository root with:

    python3 codes/run_phase_diagram.py
    python3 codes/run_spin_maps.py 4

## Reproduce

    python3 plot_FigR1.py

Needs numpy and matplotlib and the files in this folder.
