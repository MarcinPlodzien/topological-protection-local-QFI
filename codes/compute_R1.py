"""Reply figure R3: read-out scan of the submitted Fig. 2 (fermionic quantity |U_j^1+V_j^1|^2,
gamma = 1, h = 0.1J, encoding k = 1).
 R1_data.txt   : mean over tJ in [150, 200] for L = 100, 200, 400 and the zero-mode term (U_j0 U_10)^2
 R1_windows.txt: excess over the zero-mode term for L = 400 and windows [50,100], [150,200], [450,500]

Model: free transverse-field Ising chain (gamma = 1, delta = 0), h = 0.1J, J = 1, rotation R_y at site 1, all spins down.
The window mean over [T1, T2] is obtained from two [0, T] averages of free_fermion.late_time:
    mean_[T1,T2] = (T2 <F>_[0,T2] - T1 <F>_[0,T1]) / (T2 - T1),  with F = avg_y at read-out site j (j = 1..120),
which for j > 1 is the fermionic-mode expression |U_j^1 + V_j^1|^2 of the submitted Fig. 2 (no Jordan-Wigner string).
Zero-mode term: (U_j0 U_10)^2 with 0 the mode of smallest s.
R1_data.txt columns:    j avg_L100 zm_L100 avg_L200 zm_L200 avg_L400 zm_L400  (NaN for j > L)
R1_windows.txt columns: j excess_50_100 excess_150_200 excess_450_500           (L = 400)

Imports free_fermion from this folder and writes to ../data/reply/.
Run (from the repository root):  python3 codes/compute_R1.py      Runtime 61 s (2 threads).
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, '..', 'data', 'reply'); os.makedirs(OUT, exist_ok=True)
import numpy as np
from free_fermion import chiral_modes, late_time

J, gam, h, jmax = 1.0, 1.0, 0.1, 120
def window_mean(U, s, V, j, T1, T2):
    a1 = late_time(U, s, V, T1, k=0, j=j)['avg_y'] if T1 > 0 else 0.0
    a2 = late_time(U, s, V, T2, k=0, j=j)['avg_y']
    return (T2 * a2 - T1 * a1) / (T2 - T1)

cols, names = [np.arange(1, jmax + 1)], ['j']
for L in (100, 200, 400):
    U, s, V = chiral_modes(L, J, gam, h); i0 = np.argmin(s)
    n = min(jmax, L)
    avg = np.full(jmax, np.nan); zm = np.full(jmax, np.nan)
    avg[:n] = [window_mean(U, s, V, j, 150., 200.) for j in range(n)]
    zm[:n] = (U[:n, i0] * U[0, i0]) ** 2
    cols += [avg, zm]; names += [f'avg_L{L}', f'zm_L{L}']
np.savetxt(os.path.join(OUT, 'R1_data.txt'), np.column_stack(cols), header=' '.join(names), fmt='%.6e')

U, s, V = chiral_modes(400, J, gam, h); i0 = np.argmin(s); zm = (U[:jmax, i0] * U[0, i0]) ** 2
cols, names = [np.arange(1, jmax + 1)], ['j']
for T1, T2 in ((50., 100.), (150., 200.), (450., 500.)):
    cols.append(np.array([window_mean(U, s, V, j, T1, T2) for j in range(jmax)]) - zm); names.append(f'excess_{int(T1)}_{int(T2)}')
np.savetxt(os.path.join(OUT, 'R1_windows.txt'), np.column_stack(cols), header=' '.join(names), fmt='%.6e')
