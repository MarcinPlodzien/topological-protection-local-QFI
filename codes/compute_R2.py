"""Reply figure R5: edge QFI (R_y, k = j = 1, gamma = 1) near the transition in the setting of the submitted
Fig. 4: L = 50, 100, 200, 400; time-averaged QFI = mean of F(t) over the 100 times linspace(200, 300, 100) used
by the submitted generator; stationary QFI = QFI of the state averaged over tJ in [200, 300].
Output: R2_data.txt, and the edge trace at h = 0.98J in R2_trace.txt

Model: free transverse-field Ising chain (gamma = 1, delta = 0), J = 1, k = j = 1, all spins down.
F_y(t) = (sum_n U_1n^2 cos(s_n t))^2 + (sum_n V_1n U_1n sin(s_n t))^2 from free_fermion.chiral_modes.
Grid: 1 - h/J = logspace(-3, log10(0.5), 50); L = 50, 100, 200, 400.
avg = mean of F_y over the 100 times linspace(200, 300, 100) (the sampling of the submitted generator);
stat = QFI of the state averaged over [200, 300] (closed form of mean cos and mean sin).
R2_data.txt columns:  one_minus_h avg_L50 stat_L50 avg_L100 stat_L100 avg_L200 stat_L200 avg_L400 stat_L400 exact_semi_infinite
R2_trace.txt columns: tJ F_L50 F_L100 F_L200 F_L400   (h = 0.98J, tJ = linspace(0, 450, 1801))

Imports free_fermion from this folder and writes to ../data/reply/.
Run (from the repository root):  python3 codes/compute_R2.py      Runtime 4 s (2 threads).
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, '..', 'data', 'reply'); os.makedirs(OUT, exist_ok=True)
import numpy as np
from free_fermion import chiral_modes

J, gam, T1, T2 = 1.0, 1.0, 200.0, 300.0
ts = np.linspace(T1, T2, 100)
x = np.logspace(-3, np.log10(0.5), 50)                 # 1 - h/J
def mean_cos(s): return np.where(s * (T2 - T1) > 1e-12, (np.sin(s * T2) - np.sin(s * T1)) / np.maximum(s * (T2 - T1), 1e-300), np.cos(s * T1))
def mean_sin(s): return np.where(s * (T2 - T1) > 1e-12, (np.cos(s * T1) - np.cos(s * T2)) / np.maximum(s * (T2 - T1), 1e-300), np.sin(s * T1))
cols, names = [x], ['one_minus_h']
for L in (50, 100, 200, 400):
    avg, stat = [], []
    for xi in x:
        U, s, V = chiral_modes(L, J, gam, J * (1 - xi))
        aa, ba = U[0] * U[0], -V[0] * U[0]                  # G_aa(t) = sum aa cos(st), G_ba(t) = sum ba sin(st)
        F = (np.cos(np.outer(ts, s)) @ aa) ** 2 + (np.sin(np.outer(ts, s)) @ ba) ** 2
        avg.append(F.mean()); stat.append((mean_cos(s) @ aa) ** 2 + (mean_sin(s) @ ba) ** 2)
    cols += [np.array(avg), np.array(stat)]; names += [f'avg_L{L}', f'stat_L{L}']
cols.append((1 - (1 - x) ** 2) ** 2); names.append('exact_semi_infinite')
np.savetxt(os.path.join(OUT, 'R2_data.txt'), np.column_stack(cols), header=' '.join(names), fmt='%.6e')

# Edge trace F_y(t) at h = 0.98J for the four lengths, 0 <= tJ <= 450 -> R2_trace.txt
ts = np.linspace(0, 450, 1801); cols, names = [ts], ['tJ']
for L in (50, 100, 200, 400):
    U, s, V = chiral_modes(L, J, gam, 0.98 * J); aa, ba = U[0] * U[0], -V[0] * U[0]
    cols.append((np.cos(np.outer(ts, s)) @ aa) ** 2 + (np.sin(np.outer(ts, s)) @ ba) ** 2); names.append(f'F_L{L}')
np.savetxt(os.path.join(OUT, 'R2_trace.txt'), np.column_stack(cols), header=' '.join(names), fmt='%.6e')
