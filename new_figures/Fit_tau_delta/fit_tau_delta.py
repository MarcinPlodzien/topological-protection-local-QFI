"""
fit_tau_delta.py -- fit of the plateau lifetime tau(delta) of the interacting chain (reply).

Reads lifetimes.txt (copy of data/interacting/lifetimes.txt) (gamma = 1, delta > 0, L >= 16, plateau F0 > 0.1). ln(tau) is averaged over L,
with the spread over L as the uncertainty, and fitted to ln tau = a + c J/delta (prethermal form) and to
ln tau = a + p ln(1/delta) (power law), with and without the extrapolated values (tau > 5e3). For the sets with no
detected decay a lower bound on tau follows from the drop of the signal up to the end of the trace.
Output: fit_tau_delta_output.txt.   Run:  python3 fit_tau_delta.py
"""
import numpy as np, sys
import os
f=os.path.join(os.path.dirname(os.path.abspath(__file__)),"lifetimes.txt")
d=np.loadtxt(f)  # gamma h delta L F0 tau t_end drop
for h in (0.2,0.4):
    print(f"\n=== gamma=1 h={h}")
    sel=(d[:,0]==1)&(np.isclose(d[:,1],h))&(d[:,2]>0)&(d[:,3]>=16)&(d[:,4]>0.1)
    rows=[]
    for de in sorted(set(d[sel,2])):
        s=sel&np.isclose(d[:,2],de)
        tau=d[s,5]; drop=d[s,7]; tend=d[s,6]
        nd = np.all(np.abs(drop)<0.01)
        lt=np.log(tau)
        print(f"delta={de:.2f} L={d[s,3].astype(int)} tau={np.round(tau,1)} drop={np.round(drop,4)} {'n.d.' if nd else ''}")
        if nd:
            # lower bound: drop over [15,995] < max|drop| -> tau > (t_end-15)/max(drop, 1e-3)
            lb=(tend.min()-15)/max(np.abs(drop).max(),1e-3)
            rows.append((de,None,None,lb,True,False)); continue
        dag = np.all(tau>5e3)
        # weight: spread of ln(tau) over L, floored at 0.05
        rows.append((de,lt.mean(),max(lt.std(),0.05),None,False,dag))
    for incl_dag in (True,False):
        R=[r for r in rows if r[1] is not None and (incl_dag or not r[5])]
        if len(R)<3: print("too few points"); continue
        x=np.array([r[0] for r in R]); y=np.array([r[1] for r in R]); w=1/np.array([r[2] for r in R])**2
        def wfit(X):
            A=np.vstack([np.ones_like(X),X]).T; W=np.diag(w)
            cov=np.linalg.inv(A.T@W@A); p=cov@A.T@W@y
            res=y-A@p; chi2=float(res@W@res); return p,np.sqrt(np.diag(cov)),chi2,res
        pe,se,ce,re=wfit(1/x)          # ln tau = a + c/delta
        pp,sp,cp,rp=wfit(-np.log(x))   # ln tau = a + p ln(1/delta)
        print(f" fit points delta={x} (daggered extrapolated incl: {incl_dag})")
        print(f"  exp:   tau = {np.exp(pe[0]):.3g} * exp({pe[1]:.3f}+-{se[1]:.3f} J/delta)   chi2={ce:.2f} dof={len(x)-2} resid(ln)={np.round(re,2)}")
        print(f"  power: tau = {np.exp(pp[0]):.3g} * delta^-({pp[1]:.2f}+-{sp[1]:.2f})          chi2={cp:.2f} dof={len(x)-2} resid(ln)={np.round(rp,2)}")
        for r in rows:
            if r[4]:
                de=r[0]; print(f"   n.d. delta={de}: lower bound tau>{r[3]:.2g};  exp-fit predicts {np.exp(pe[0]+pe[1]/de):.2g}, power predicts {np.exp(pp[0]-pp[1]*np.log(de)):.2g}")

# Unweighted fits (no uncertainty model): slope and rms residual of ln(tau), for the same points
print("\n=== unweighted fits of ln(tau), L >= 16, F0 > 0.1, decaying sets")
for h in (0.2, 0.4):
    sel = (d[:, 0] == 1) & np.isclose(d[:, 1], h) & (d[:, 2] > 0) & (d[:, 3] >= 16) & (d[:, 4] > 0.1)
    xs, ys = [], []
    for de in sorted(set(d[sel, 2])):
        s = sel & np.isclose(d[:, 2], de)
        if np.all(np.abs(d[s, 7]) < 0.01): continue
        xs.append(de); ys.append(np.log(d[s, 5]).mean())
    x, y = np.array(xs), np.array(ys)
    for name, X in (('exp  c (J/delta)', 1 / x), ('power p (1/delta)', -np.log(x))):
        A = np.vstack([np.ones_like(X), X]).T
        p_, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ p_
        print(f"  h={h}: {name}: slope {p_[1]:.2f}, rms residual of ln(tau) {np.sqrt(r @ r / len(r)):.3f}")

# Unweighted fits using only resolved lifetimes (decay detected and tau <= 5e3, i.e. not extrapolated)
print("\n=== unweighted fits, resolved lifetimes only (no n.d., no tau > 5e3)")
for h in (0.2, 0.4):
    sel = (d[:, 0] == 1) & np.isclose(d[:, 1], h) & (d[:, 2] > 0) & (d[:, 3] >= 16) & (d[:, 4] > 0.1)
    xs, ys = [], []
    for de in sorted(set(d[sel, 2])):
        s = sel & np.isclose(d[:, 2], de)
        if np.all(np.abs(d[s, 7]) < 0.01) or np.all(d[s, 5] > 5e3): continue
        xs.append(de); ys.append(np.log(d[s, 5]).mean())
    x, y = np.array(xs), np.array(ys)
    print(f"  h={h}: deltas {x}")
    for name, X in (('exp  c (J/delta)', 1 / x), ('power p (1/delta)', -np.log(x))):
        A = np.vstack([np.ones_like(X), X]).T
        p_, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ p_
        print(f"    {name}: slope {p_[1]:.2f}, rms residual of ln(tau) {np.sqrt(r @ r / len(r)):.3f}")
