#!/usr/bin/env python3
"""
plot_long_traces_logt.py -- the diagnostic figure of plot_long_traces.py with a logarithmic time axis
(not a figure of the paper).

Top row: F_y(t), running average over 10/J, log-log, for L = 16, 18, 20 (traces to tJ = 10^4, run_long_traces.py)
and L = 22 (to tJ = 10^3, run_interacting_long.py), at (h, delta) = (0.2J, 0.2J) and (0.4J, 0.1J), gamma = 1.
Dashed vertical lines: Heisenberg-time estimates of data/interacting/heisenberg_time.txt inside the window.
Bottom row: local decay time tau_loc(t) = -1 / (d ln F_y / dt) from linear fits of ln F_y on windows [t, 1.5 t]
(log-spaced t), which suits the logarithmic axis.
Output: figures/long_traces_check_logt.pdf and .png.   Run: python3 codes/plot_long_traces_logt.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'data')
COL = {16: '#2a78d6', 18: '#eb6834', 20: '#1baf7a', 22: '#eda100'}   # fixed categorical order, one colour per L
SETS = ((0.2, 0.2), (0.4, 0.1))

def smooth(t, y, w=10.0):
    n = max(1, int(round(w / (t[1] - t[0])))); return t[n - 1:], np.convolve(y, np.ones(n) / n, 'valid')

def trace(L, h, d):
    f = (os.path.join(D, 'interacting_long', f'trace_L{L}_g1.00_h{h:.2f}_d{d:.2f}_t10000.txt') if L < 22 else
         os.path.join(D, 'interacting', f'trace_L{L}_g1.00_h{h:.2f}_d{d:.2f}.txt'))
    return np.loadtxt(f)

ht = np.loadtxt(os.path.join(D, 'interacting', 'heisenberg_time.txt'))
fig, ax = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
for c, (h, d) in enumerate(SETS):
    for L in (16, 18, 20, 22):
        x = trace(L, h, d); t, y = smooth(x[:, 0], x[:, 1])
        ax[0, c].loglog(t, y, color=COL[L], lw=1.4, label=f'L = {L}' + (' ($tJ\\leq10^3$)' if L == 22 else ''))
        th = ht[(abs(ht[:, 0] - h) < 1e-9) & (abs(ht[:, 1] - d) < 1e-9) & (ht[:, 2] == L), 8]
        if len(th) and th[0] < t[-1]:
            ax[0, c].axvline(th[0], color=COL[L], ls='--', lw=1.0)
            ax[0, c].text(th[0], 0.5, ' $t_H$', color='#52514e', fontsize=8, va='bottom')
        tc, tl = [], []
        for a in np.geomspace(15, t[-1] / 1.5, 30):
            k = (t >= a) & (t <= 1.5 * a)
            s = np.polyfit(t[k], np.log(y[k]), 1)[0]
            tc.append(1.22 * a); tl.append(-1 / s if s < 0 else np.nan)
        ax[1, c].loglog(tc, tl, 'o-', color=COL[L], ms=3, lw=1.0)
    ax[0, c].set_title(f'$h = {h}J$, $\\delta = {d}J$', fontsize=10)
    ax[0, c].set_ylim(1e-3, 1.0); ax[1, c].set_ylim(1e2, 1e5); ax[1, c].set_xlim(8, 1.2e4)
    ax[1, c].set_xlabel('$tJ$')
    for a in ax[:, c]:
        a.grid(True, which='major', color='#e6e5e0', lw=0.6); a.set_axisbelow(True)
        for s in ('top', 'right'):
            a.spines[s].set_visible(False)
ax[0, 0].set_ylabel('$F_y(t)$ (running average over $10/J$)')
ax[1, 0].set_ylabel(r'local decay time $\tau_{\rm loc}J$ (windows $[t, 1.5t]$)')
ax[0, 0].legend(frameon=False, fontsize=8, loc='lower left')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, '..', 'figures', f'long_traces_check_logt.{ext}'), dpi=150)
print('written figures/long_traces_check_logt.{pdf,png}')
