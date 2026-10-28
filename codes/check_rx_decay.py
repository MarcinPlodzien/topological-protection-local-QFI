"""Decay of the edge propagators at L=1500: |G| t^{3/2} and the log-log slope of the envelope maxima.
G^{aa}_{11} - N_1^2 enters R_y; G^{bb}_{11}, G^{ab}_{11} and I_x = (G^{ab})^2 + (G^{bb})^2 enter R_x.
Run: python3 codes/check_rx_decay.py > data/free/rx_decay_L1500.txt

Model: free XY chain (delta = 0), J = 1, L = 1500 (no reflection from the far end for tJ <= 640), k = j = 1;
(gamma, h) = (1, 0.5), (1, 0.2), (0.5, 0.5), (0.5, 1.5).  Propagators on the grid tJ = linspace(10, 640, 126001):
G^aa_11 - N_1^2 (zero-mode part removed when the smallest s < 1e-6), G^bb_11, G^ab_11, and I_x = (G^ab)^2 + (G^bb)^2.
For each: range of the local maxima of |G| t^p (p = 3/2, and p = 3 for I_x) on [10,40], [40,160], [160,640], and the
log-log slope of the local maxima of |G| between tJ = 40 and 640.  Quoted in the reply,
and in the reply (axis asymmetry).
Memory: the two 126001 x 1500 cos/sin matrices take 1.5 GB each; peak 6.1 GB.  Runtime 40 s (2 threads).
"""
import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from free_fermion import chiral_modes
from scipy.signal import argrelmax
L=1500
t=np.linspace(10,640,126001)
def env(y):
    i=argrelmax(y)[0]; return t[i], y[i]
for gam,h in ((1.0,0.5),(1.0,0.2),(0.5,0.5),(0.5,1.5)):
    U,s,V=chiral_modes(L,1.0,gam,h)
    z=np.argmin(s); print(f"gamma={gam} h={h}: smallest s={s[z]:.2e}, U0={U[0,z]:.3f}, V0={V[0,z]:.2e}")
    wbb=V[0]**2; wab=U[0]*V[0]; waa=U[0]**2
    C=np.cos(np.outer(t,s)); S=np.sin(np.outer(t,s))
    Gbb=C@wbb; Gab=S@wab; Gaa=C@waa-(waa[z] if s[z]<1e-6 else 0)
    Ix=Gab**2+Gbb**2
    for name,G,p in (('Gaa-N1^2',Gaa,1.5),('Gbb',Gbb,1.5),('Gab',Gab,1.5),('I_x',Ix,3.0)):
        tt,yy=env(np.abs(G)*t**p)
        out=[]
        for a,b in ((10,40),(40,160),(160,640)):
            m=(tt>=a)&(tt<=b); out.append(f"[{a},{b}] {yy[m].min():.3f}-{yy[m].max():.3f}")
        # log-log slope of the envelope maxima of |G| between 40 and 640
        tt2,yy2=env(np.abs(G)); m=(tt2>=40)&(tt2<=640)
        sl=np.polyfit(np.log(tt2[m]),np.log(yy2[m]),1)[0]
        print(f"  {name:9s} |G| t^{p}: "+"; ".join(out)+f"; slope of envelope of |G| (40..640) = {sl:.3f}")
