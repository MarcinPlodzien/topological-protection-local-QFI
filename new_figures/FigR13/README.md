# FigR13: Interacting chain with a trivial zero-energy edge state

FigR13 of the reply to the referees (answer to Referee B, point 1).

## Panels

(a) Edge QFI versus time, running average over 10/J (time at the window centre), L = 16, γ = 1, δ = 0.1J: trivial chain with h_1 = 0 and h_j = 1.8J for j ≥ 2 (red) and topological chain with h = 0.2J (blue); rotation about y (solid) and about x (dashed).

(b) Stationary QFI, the QFI of the state averaged over [0, T], versus T, full diagonalization, L = 12, for the same two chains; grey: the trivial chain at δ = 0.

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| trivial_L{14,16}_h1{0.00,0.30}_d{...}.txt | columns tJ, F_y, F_x, output every 0.5/J to tJ = 1000; h_1 = 0 or 0.3J, δ = 0 … 0.2J | run_trivial_interacting.py |
| trace_L16_g1.00_h0.20_d0.10.txt | topological chain, columns tJ, F_y, F_x (copy of ../../data/interacting/) | ../../codes/run_interacting.py |
| stationary_vs_T_<case>.txt | columns T, I_st(R_y), I_st(R_x), 21 values of T between 1/J and 10^4/J | stationary_vs_T.py |
| stationary_L12.txt | stationary QFI at TJ = 10^2, 10^3, 10^4 for the topological and trivial chains, δ = 0, 0.1J, 0.2J | python3 stationary_qfi_ed.py > stationary_L12.txt |
| stationary_qfi_ed.py | QFI of the time-averaged state of the edge spin by full diagonalization, with the window kernel | |

## Reproduce

    python3 run_trivial_interacting.py 8      # traces, 8 processes
    python3 stationary_vs_T.py                # panel (b) data, L = 12
    python3 stationary_qfi_ed.py > stationary_L12.txt
    python3 plot_FigR13.py                    # writes FigR13.pdf

The plot script needs numpy and matplotlib and the files in this folder.
