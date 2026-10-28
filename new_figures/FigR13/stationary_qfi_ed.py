"""Stationary QFI (QFI of the time-averaged edge state over [0,T]) by full diagonalization, exact time average.
H = -(J/2) sum[(1+g)XX+(1-g)YY] - sum h_j Z + delta sum ZZ, all down, rotation of edge spin about y or x."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "8")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "codes"))
import ed_interacting as ed
def sigma_edge(L):
    N=1<<L; idx=np.arange(N); up=(idx&1)==1          # bit 0 = edge spin, 1 = up
    flip=idx^1
    X=np.zeros((N,N)); X[flip,idx]=1.0
    Y=np.zeros((N,N),complex); Y[flip,idx]=np.where(up,-1j,1j)   # sigma^y|dn>=i|up>, |up>-> -i|dn>
    Z=np.diag(np.where(up,1.0,-1.0))
    return X,Y,Z
def kern(w,T):  # (1/T) int_0^T exp(i w t) dt
    x=w*T; out=np.ones_like(x,dtype=complex); m=np.abs(x)>1e-12
    out[m]=(np.exp(1j*x[m])-1)/(1j*x[m]); return out
def stationary(L,g,h,delta,Ts):
    diag,off=ed.build_H(L,1.0,g,h,delta); E,V=np.linalg.eigh(off.toarray()+np.diag(diag))
    psi0,phi=ed.initial_states(L,0); a=V.T@psi0; b=V.T@phi
    S=[V.T@M@V for M in sigma_edge(L)]
    W=E[:,None]-E[None,:]                    # <psi(t)|O|chi(t)> = sum_mn a_m* O_mn chi_n e^{i(E_m-E_n)t}
    res={}
    for T in Ts:
        K=kern(W,T)
        m=np.array([np.real(np.conj(a)@(Om*K)@a) for Om in S])
        c=np.array([np.conj(a)@(Om*K)@b for Om in S])
        out=[]
        for pref in (-0.5,-0.5j):
            dm=2.0*np.real(pref*c)
            out.append(dm@dm+(m@dm)**2/max(1-m@m,1e-14))
        res[T]=out
    # instantaneous check at t=0+: QFI(t=0) must be 1 for R_y
    return res
if __name__=="__main__":
    L=12; Ts=(1e2,1e3,1e4)
    cases=[("topological h=0.2",lambda: np.full(L,0.2)),("topological h=0.4",lambda: np.full(L,0.4)),
           ("trivial h=1.8, h1=0",lambda: np.r_[0.0,np.full(L-1,1.8)]),("trivial h=1.8, h1=0.3",lambda: np.r_[0.3,np.full(L-1,1.8)])]
    print("# L=12, gamma=1; stationary QFI R_y / R_x for TJ =",Ts)
    for name,hf in cases:
        for d in (0.0,0.1,0.2):
            r=stationary(L,1.0,hf(),d,Ts)
            print(f"{name:22s} delta={d:4.2f}  "+"   ".join(f"T={T:.0e}: {r[T][0]:.4f}/{r[T][1]:.4f}" for T in Ts),flush=True)
