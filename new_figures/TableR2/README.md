# TableR2: Plateau value and lifetime of the edge QFI in the interacting chain

TableR2 of the reply to the referees.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| trace_L{L}_g{γ}_h{h}_d{δ}.txt, L = 12…18 (56 files) | columns tJ, F_y, F_x; the traces of run_interacting.py for these sizes without the four at h = 0.4J, δ = 0.4J (F0 ≈ 0.02, not in the table) | ../../codes/run_interacting.py |
| trace_L20_g1.00_h{h}_d{δ}.txt (11 files) | same; the twelve γ = 1 sets of the table without h = 0.4J, δ = 0.4J | ../../codes/run_interacting_long.py |
| trace_L22_g1.00_h{h}_d{δ}.txt (8 files) | same; h = 0.2J with δ = 0.15, 0.2, 0.3 and h = 0.4J with δ = 0.05, 0.1, 0.15, 0.2, 0.3 | ../../codes/run_interacting_long.py |
| lifetimes.txt | columns γ, h, δ, L, F0, τ, t_end, drop for all 85 traces of ../../data/interacting/ (including the h = 0 XYZ traces and h = 0.4J, δ = 0.4J); written with 4 decimals of F0, 5 significant digits of τ, 1 decimal of t_end and 4 decimals of drop | ../../codes/make_figures.py (function lifetimes_table, called by fig4) |

## Data

Generated from the repository root with:

    python3 codes/make_tableR2.py

## Reproduce

    python3 plot_TableR2.py

Needs numpy and matplotlib and the files in this folder.
