#!/usr/bin/env python3
"""
plot_TableR1.py -- prints TableR1 of the reply from the files cases_*.txt in this folder.

Edge response for R_y / R_x at gamma = 1, L = 200, k = j = 1. Columns:
  eps_low/J            lowest single-particle energy (smallest singular value of B), read from the file header
  avg QFI (TJ = 1e4)   time-averaged QFI over the window [0, T]
  I_st (TJ = 31.6)     stationary QFI (QFI of the time-averaged state); the grid point is T = 10^1.5 = 31.62
  I_st (TJ = 1e6)      the same at T = 1e6
The averaging times are grid points of the file (np.logspace(0.5, 10, 77)); no interpolation is done.

Run:  python3 plot_TableR1.py      (prints the table; well under 1 s).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = (('topological', 'Topological, h=0.5J'), ('trivial_clean', 'Trivial, h=1.8J'), ('trivial_impurity', 'Trivial, h_1=0.3J'),
        ('trivial_finetuned', 'Trivial, h_1=0'), ('island_n8', 'Island, n=8'), ('island_n16', 'Island, n=16'),
        ('smooth_ramp', 'Ramp over 30 sites'))


def at(c, T):
    """row of the table c whose averaging time (column 0) is the grid point closest to T on a log scale"""
    return c[np.argmin(np.abs(np.log10(c[:, 0]) - np.log10(T)))]


if __name__ == "__main__":
    print(f"{'Configuration':22s} {'eps_low/J':>10s}   {'avg I (TJ=1e4)':>17s}   {'I_st (TJ=31.6)':>17s}   {'I_st (TJ=1e6)':>21s}")
    print("                                       R_y / R_x           R_y / R_x               R_y / R_x")
    for nm, lab in ROWS:
        f = os.path.join(HERE, f'cases_{nm}.txt')
        with open(f) as fh:
            eps = float(fh.readline().split(':')[1].split()[0])      # '# L=200 gamma=1.0; three lowest |omega|: e0 e1 e2'
        c = np.loadtxt(f)                                            # columns: T stat_y stat_x avg_y avg_x
        r4, r32, r6 = at(c, 1e4), at(c, 10 ** 1.5), at(c, 1e6)
        print(f"{lab:22s} {eps:10.2e}   {r4[3]:7.3f} / {r4[4]:7.3f}   {r32[1]:7.3f} / {r32[2]:7.3f}   {r6[1]:9.2e} / {r6[2]:9.2e}")
