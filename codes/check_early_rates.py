#!/usr/bin/env python3
"""
check_early_rates.py -- decay rate of the edge QFI before the first return of excitations from the far end
(reply).

For each gamma = 1 trace in data/interacting/ with a resolved decay (h = 0.2J, delta = 0.3J, 0.4J; h = 0.4J,
delta = 0.2J, 0.3J), the running average of F_y over 10/J (as in make_figures.smooth) is fitted by ln F_y = a - t/tau on
15 <= tJ <= 2L/v with the free-fermion velocity v = 2 min(h, J). Columns: h, delta, L, tau_early, end of the window,
tau of TableR2 (lifetimes.txt).
Run: python3 codes/check_early_rates.py > data/interacting/early_rates.txt   (seconds)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data', 'interacting')
src = open(os.path.join(HERE, 'make_figures.py')).read(); ns = {'np': np}
i = src.index('def smooth'); exec(src[i:src.index('\ndef ', i + 5)], ns); smooth = ns['smooth']
T = np.loadtxt(os.path.join(D, 'lifetimes.txt'))
print("# gamma=1; tau from a fit of ln F_y on 15 <= tJ <= 2L/v, v = 2 min(h,J); columns: h delta L tau_early t_window tau_TableR2")
for h, d in ((0.2, 0.3), (0.2, 0.4), (0.4, 0.2), (0.4, 0.3)):
    for L in (16, 18, 20, 22):
        f = os.path.join(D, f'trace_L{L}_g1.00_h{h:.2f}_d{d:.2f}.txt')
        if not os.path.exists(f):
            continue
        x = np.loadtxt(f); t, y = smooth(x[:, 0], x[:, 1], 10.0)
        tw = 2 * L / (2 * min(h, 1.0)); k = (t >= 15) & (t <= tw)
        tau = -1 / np.polyfit(t[k], np.log(y[k]), 1)[0]
        m = (T[:, 0] == 1) & (abs(T[:, 1] - h) < 1e-9) & (abs(T[:, 2] - d) < 1e-9) & (T[:, 3] == L)
        print(f"{h:.2f} {d:.2f} {L} {tau:.1f} {tw:.0f} {T[m, 5][0]:.1f}")
