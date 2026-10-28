# FigR9: Interacting chain: edge QFI traces and lifetime of the plateau

FigR9 of the reply to the referees.

## Panels

(a) Running average over 10/J of F_y(t) (solid) and F_x(t) (dotted), h = 0.2J, L = 18, δ/J = 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, for 5 ≤ tJ ≤ 1000 (log scale). Open circles: L = 20 for all six δ, 14 points per curve on a geometric grid.

(b) Lifetime τ of the F_y plateau versus δ, h = 0.2J (circles) and 0.4J (squares), L = 16 (large open symbols), 18 (small filled), 20 (stars), 22 (crosses). Rows with plateau value F0 ≤ 0.1 or δ ≤ 0 are not drawn. Points are drawn when the running average dropped by more than 1 % over the trace; otherwise an arrow marks "no decay detected" (for L = 18, or for L = 16 when no L = 18 trace exists). Grey band: τ > 5 × 10³/J. All drawn points in the band have a drop of at most 7.0 % up to tJ = 10³ (largest: 0.0696), all points below it at least 27 % (smallest: 0.274), so τ in the band is an extrapolation of the fitted slope. The values near 2 × 10³/J at (h, δ) = (0.2J, 0.2J) and (0.4J, 0.1J) are tested by the traces to tJ = 10⁴ of FigR11.

Lifetime definition (function lifetime(), copied unchanged from ../../codes/make_figures.py): running average of F_y over 10/J; F0 = mean of the running average over 5 ≤ tJ ≤ 15; t_end = first time after tJ = 15 at which the running average is below 0.15 F0, or the end of the trace; τ = -1/slope of a linear least-squares fit of the logarithm of the running average on 15 ≤ tJ ≤ t_end; drop = 1 - (last value of the running average)/F0.

## Files in this folder

trace_L{L}_g1.00_h{h}_d{δ}.txt, 41 files: L = 16, 18 and 20 for h = 0.2J with δ = 0.05, 0.1, 0.15, 0.2, 0.3, 0.4 and for h = 0.4J with δ = 0.05, 0.1, 0.15, 0.2, 0.3; L = 22 for h = 0.2J with δ = 0.15, 0.2, 0.3 and for h = 0.4J with δ = 0.05, 0.1, 0.15, 0.2, 0.3. Columns: tJ, F_y, F_x; tJ = 0.5, 1.0, …, 1000. The header records the final norm error of the propagation (at most 1.0 × 10⁻¹² in modulus over all 85 traces of ../../data/interacting/).

The traces at h = 0.4J, δ = 0.4J (L = 12…20) exist in ../../data/interacting/ but have F0 ≈ 0.02 and are excluded from panel (b) by the condition F0 > 0.1; they are not copied here. The plot script computes the lifetime table from the traces in this folder only; ../../codes/make_figures.py uses all traces in ../../data/interacting/ and writes ../../data/interacting/lifetimes.txt. The rows that enter the figure are the same.

## Data

Generated from the repository root with:

    python3 codes/run_interacting.py 12
    python3 codes/run_interacting_long.py 20 2 6
    python3 codes/run_interacting_long.py 22 2 6 8

## Reproduce

    python3 plot_FigR9.py

Needs numpy and matplotlib and the files in this folder.
