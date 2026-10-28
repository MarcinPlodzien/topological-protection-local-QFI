# FigR2: Single-spin QFI with the Jordan-Wigner string against the fermionic-mode expression

FigR2 of the reply to the referees.

## Panels

(a) Profiles over j at tJ = 30 (the grid point nearest to 30) for the topological chain, γ = 0.3, h = 0.5J: single-spin QFI (circles) and |U_j^1 + V_j^1|² = G[a_j,a_1]² + G[b_j,a_1]² (squares). Values below 10⁻¹² are clipped to 10⁻¹² for the log scale.

(b) The same for the trivial chain, γ = 0.8, h = 1.8J.

(c) Both quantities summed over the 60 stored sites versus time; solid w = 1 (topological), dashed w = 0 (trivial).

## Files in this folder

| File | Content | Produced by |
|---|---|---|
| spinmap_times.txt | the 121 times tJ | ../../codes/run_spin_maps.py |
| spinmap_topological_Ry.txt | 121 × 60, rows time, columns j = 1…60: exact single-spin QFI, R_y, γ = 0.3, h = 0.5J | ../../codes/run_spin_maps.py |
| spinmap_topological_fermionic.txt | 121 × 60: \|U_j^1 + V_j^1\|², same parameters | ../../codes/run_spin_maps.py |
| spinmap_trivial_Ry.txt | as above for γ = 0.8, h = 1.8J | ../../codes/run_spin_maps.py |
| spinmap_trivial_fermionic.txt | as above for γ = 0.8, h = 1.8J | ../../codes/run_spin_maps.py |

## Data

Generated from the repository root with:

    python3 codes/run_spin_maps.py [nproc]
    python3 codes/run_spin_maps.py 2

## Reproduce

    python3 plot_FigR2.py

Needs numpy and matplotlib and the files in this folder.
