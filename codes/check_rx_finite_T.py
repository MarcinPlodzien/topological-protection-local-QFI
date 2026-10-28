"""Finite-T correction to the stationary R_x QFI in the topological chain (reply, regime (ii)).

Reads data/free/cases_topological.txt (gamma = 1, h = 0.5J, L = 200; written by run_free.py) and prints
I_st,(x) * (JT)^2 versus T. A constant shows that I_st,(x) = c/(JT)^2 + O(e^{-4L/xi}).
Run: python3 codes/check_rx_finite_T.py > data/free/rx_finite_T.txt   (seconds)"""
import os
import numpy as np

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'free')
x = np.loadtxt(os.path.join(D, 'cases_topological.txt'))
print("# from cases_topological.txt (gamma=1, h=0.5J, L=200); columns: TJ, I_st_x, I_st_x*(TJ)^2")
for T, sx in zip(x[:, 0], x[:, 2]):
    if 10 <= T <= 750:
        print(f"{T:.4g} {sx:.4e} {sx * T * T:.4f}")
v = x[(x[:, 0] >= 10) & (x[:, 0] <= 750)]
c = v[:, 2] * v[:, 0] ** 2
print(f"# 10 <= TJ <= 750: mean {c.mean():.4f}, min {c.min():.4f}, max {c.max():.4f}")
