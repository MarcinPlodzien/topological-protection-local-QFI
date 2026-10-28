# FigR11: Long traces of the interacting chain to tJ = 10^4

FigR11 of the reply to the referees.

## Panels

(a) Running average of F_y(t), log-log, (h, δ) = (0.2J, 0.2J): L = 16 (dark grey), 18 (black), 20 (magenta) to tJ = 10⁴ and L = 22 (indigo) to tJ = 10³. Axis range 10⁻³ ≤ F_y ≤ 1.2.

(b) The same for (h, δ) = (0.4J, 0.1J).

(c, d) Local decay time τ_loc(t) = −[d ln F_y/dt]⁻¹ for the curves of (a) and (b): for 30 values of t0 log-spaced between 15/J and t_last/1.5, a linear fit of the logarithm of the running average on [t0, 1.5 t0], plotted at 1.22 t0; windows with a positive slope are not drawn. Axis range 10² ≤ τ_loc J ≤ 10⁵, 8 ≤ tJ ≤ 1.2 × 10⁴. A horizontal line corresponds to an exponential decay.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| trace_L{16,18,20}_g1.00_h0.20_d0.20_t10000.txt, trace_L{16,18,20}_g1.00_h0.40_d0.10_t10000.txt | columns tJ, F_y, F_x; tJ = 0.5, 1.0, …, 10⁴ (20000 rows); the header records the final norm error (7.5 × 10⁻¹⁴ to 6.8 × 10⁻¹² in modulus) | ../../codes/run_long_traces.py |
| trace_L22_g1.00_h0.20_d0.20.txt, trace_L22_g1.00_h0.40_d0.10.txt | columns tJ, F_y, F_x; to tJ = 1000 | ../../codes/run_interacting_long.py |

The plot script uses a *_t10000.txt file where it exists and the tJ ≤ 1000 trace otherwise, as figR11() in ../../codes/make_figures.py. Over tJ ≤ 1000 each tJ = 10⁴ trace agrees with the tJ ≤ 1000 trace of the same size in ../../data/interacting/ to 1.0 × 10⁻⁹ or better (the last stored digit).

## Data

Generated from the repository root with:

    python3 codes/run_long_traces.py 2 6
    python3 codes/run_interacting_long.py 22 2 6 8

## Reproduce

    python3 plot_FigR11.py

Needs numpy and matplotlib and the files in this folder.
