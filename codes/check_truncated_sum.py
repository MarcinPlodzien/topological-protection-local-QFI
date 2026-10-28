"""Finite-length contribution near the disorder window of the clean-trivial chain (reply, FigR7(d)).

gamma = 1, h = 1.5J, fields h_j uniform in [h - Delta, h + Delta]. For each realization the stationary edge QFI of the
semi-infinite chain, I_st = [sum_{j>=1} prod_{l<j} (h_l/J)^2]^(-2), is evaluated with the sum cut at j <= L = 200,
which mimics a chain of 200 sites. Output: Delta, mean, standard error of the mean, median, number of realizations.
Run: python3 codes/check_truncated_sum.py > data/disorder/truncated_sum_L200_h1.50.txt   (about 1 min, 1 GB)"""
import numpy as np

L, h, N, CHUNK = 200, 1.5, 400000, 20000
print(f"# gamma=1 h={h} sum cut at L={L}; {N} realizations per Delta; seed 20261015 + 1000*Delta")
print("# Delta mean sem median N")
for D in (1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.10, 2.25):
    rng = np.random.default_rng(20261015 + int(round(1000 * D)))
    vals = []
    for _ in range(N // CHUNK):
        hs = rng.uniform(h - D, h + D, (CHUNK, L))
        lr = 2.0 * np.log(np.abs(hs) + 1e-300)
        c = np.concatenate((np.zeros((CHUNK, 1)), np.cumsum(lr, axis=1)[:, :-1]), axis=1)
        m = c.max(axis=1, keepdims=True)
        vals.append(np.exp(-2.0 * (m[:, 0] + np.log(np.exp(c - m).sum(axis=1)))))
    v = np.concatenate(vals)
    print(f"{D:.2f} {v.mean():.5f} {v.std(ddof=1)/np.sqrt(v.size):.5f} {np.median(v):.3e} {v.size}", flush=True)
