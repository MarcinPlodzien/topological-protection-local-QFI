#!/usr/bin/env python3
"""
check_efold_times.py -- time at which the edge QFI has fallen to 1/e of its plateau, read directly from the traces to
tJ = 10^4 (no fit; reply).

F0 = mean of F_y over 5 <= tJ <= 15 (raw samples); running average of F_y over the preceding 10/J; the e-folding time
is the first tJ at which the running average falls below F0/e. The smoothed signal fluctuates, so it crosses F0/e
several times within about 10-20/J; the first crossing is reported.
Data: data/interacting_long/trace_L{16,18,20}_g1.00_h*_d*_t10000.txt (run_long_traces.py).
Columns: h, delta, L, F0, t(F0/e), tau of TableR2 (lifetimes.txt, fit of ln F_y on 15 <= tJ <= t_end), mean F_y and
mean F_x over 9000 <= tJ <= 10000, smallest and largest local decay time on the windows [t, 1.5t] of FigR11 for
100 <= t <= 6667 (same running average as make_figures.smooth).
Run: python3 codes/check_efold_times.py > data/interacting_long/efold_times.txt   (seconds)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data')
T = np.loadtxt(os.path.join(D, 'interacting', 'lifetimes.txt'))
print("# columns: h delta L F0 t_(F0/e) tau_TableR2 Fy_late Fx_late tauloc_min tauloc_max")
for h, d in ((0.2, 0.2), (0.4, 0.1)):
    for L in (16, 18, 20):
        x = np.loadtxt(os.path.join(D, 'interacting_long', f'trace_L{L}_g1.00_h{h:.2f}_d{d:.2f}_t10000.txt'))
        t, y = x[:, 0], x[:, 1]; dt = t[1] - t[0]; n = int(round(10.0 / dt))
        ra = np.convolve(y, np.ones(n) / n, 'valid'); tr = t[n - 1:]          # average over the preceding 10/J
        F0 = y[(t >= 5) & (t <= 15)].mean()
        te = tr[np.argmax(ra < F0 / np.e)]
        m = (T[:, 0] == 1) & (abs(T[:, 1] - h) < 1e-9) & (abs(T[:, 2] - d) < 1e-9) & (T[:, 3] == L)
        late = x[:, 0] >= 9000
        tl = []
        for t0 in np.geomspace(15, tr[-1] / 1.5, 30):          # windows of make_figures.figR11
            k = (tr >= t0) & (tr <= 1.5 * t0)
            if t0 >= 100:
                sl = np.polyfit(tr[k], np.log(ra[k]), 1)[0]
                if sl < 0:
                    tl.append(-1.0 / sl)
        print(f"{h:.2f} {d:.2f} {L} {F0:.4f} {te:.1f} {T[m, 5][0]:.1f} {x[late, 1].mean():.4f} {x[late, 2].mean():.4f} "
              f"{min(tl):.0f} {max(tl):.0f}")
