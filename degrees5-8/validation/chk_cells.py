"""Numerical tests of the per-cell inequalities used by verify_flint.
Upper-index side: finite-element condensed energy of Q_{gL} minus the
verifier's chord cell matrix. Lower-index side: the verifier's lower cell
matrix minus the energy of the piecewise-linear trial function under Q_{gU}.
Finite-element condensation and quadrature do not provide rigorous enclosures.
"""
import json,gzip,sys,importlib.util
from pathlib import Path
import mpmath as mp
mp.mp.dps=50
root=str(Path(__file__).resolve().parents[1])+'/'
spec=importlib.util.spec_from_file_location('vf',root+'verify_flint.py');vf=importlib.util.module_from_spec(spec);spec.loader.exec_module(vf)
def load(n):
    r=json.loads(open(root+f'profiles/profile_n{n}.json').read());hs=vf.profiles(r);return r,hs
def mpq(q):return mp.mpf(int(q.numerator))/mp.mpf(int(q.denominator))
def hpoly(H):
    co=[mpq(c) for c in H.coeffs()]
    def ev(y):
        s=mp.mpf(0)
        for c in reversed(co):s=s*y+c
        return s
    return ev
# 4-pt Gauss-Legendre on [0,1]
def gl4():
    import numpy as np
    x,w=np.polynomial.legendre.leggauss(4)
    # refine nodes to 50 digits with Newton on P4
    nodes=[]
    for x0 in x:
        t=mp.mpf(float(x0))
        for _ in range(8):
            t=t-mp.legendre(4,t)/ (4*(t*mp.legendre(4,t)-mp.legendre(3,t))/(t*t-1))
        nodes.append(t)
    wts=[2*(1-t*t)/(4*mp.legendre(3,t))**2 for t in nodes]
    return [(t+1)/2 for t in nodes],[w/2 for w in wts]
GX,GW=gl4()
def Bmat(n,Hev,k,j,t):
    x=mp.e**(2*t);y=x/(n*n+x);g2=y**n*Hev(y)**2
    return mp.matrix([[k*k-x+2*x*g2,x*g2],[x*g2,j*j-x+2*x*g2]])
def fe_condense(n,Hev,k,j,ta,ell,N):
    h=ell/N;dim=2*(N+1);K=mp.zeros(dim,dim)
    for e in range(N):
        t0=ta+e*h
        for xg,wg in zip(GX,GW):
            t=t0+xg*h;B=Bmat(n,Hev,k,j,t);phi=[1-xg,xg]
            for p in range(2):
                for q in range(2):
                    for c in range(2):
                        for d in range(2):
                            K[2*(e+p)+c,2*(e+q)+d]+=wg*h*phi[p]*phi[q]*B[c,d]
        for p in range(2):
            for q in range(2):
                s=(1 if p==q else -1)/h
                for c in range(2):K[2*(e+p)+c,2*(e+q)+c]+=s
    bd=[0,1,dim-2,dim-1];it=[i for i in range(dim) if i not in bd]
    Kbb=mp.matrix(4,4);Kbi=mp.matrix(4,len(it));Kii=mp.matrix(len(it),len(it))
    for a,ia in enumerate(bd):
        for b,ib in enumerate(bd):Kbb[a,b]=K[ia,ib]
        for b,ib in enumerate(it):Kbi[a,b]=K[ia,ib]
    for a,ia in enumerate(it):
        for b,ib in enumerate(it):Kii[a,b]=K[ia,ib]
    cols=[mp.lu_solve(Kii,Kbi.T[:,a]) for a in range(4)]
    X=mp.matrix(len(it),4)
    for a in range(4):
        for r_ in range(len(it)):X[r_,a]=cols[a][r_]
    return Kbb-Kbi*X
def lin_energy(n,Hev,k,j,ta,ell):
    E=mp.zeros(4,4)
    for xg,wg in zip(GX,GW):
        t=ta+xg*ell;B=Bmat(n,Hev,k,j,t);phi=[1-xg,xg]
        for p in range(2):
            for q in range(2):
                for c in range(2):
                    for d in range(2):E[2*p+c,2*q+d]+=wg*ell*phi[p]*phi[q]*B[c,d]
    for p in range(2):
        for q in range(2):
            for c in range(2):E[2*p+c,2*q+c]+=(1 if p==q else -1)/ell
    return E
def mineig(A):
    A=(A+A.T)/2;return min(mp.eigsy(A,eigvals_only=True))
def run(n,ms,cells,N=24):
    r,hs=load(n);HL=hpoly(hs['lower']);HU=hpoly(hs['upper'])
    c=json.loads(gzip.decompress(open(root+f'certificates/upper_n{n}_h{"0.001" if n==7 else "0.002"}.json.gz','rb').read()))
    xs=[vf.Q(s) for s in c['x_grid']]
    for m in ms:
        k=abs(n-m);j=n+m
        for i in cells:
            a,b=xs[i],xs[i+1];ta=mp.log(mpq(a))/2;tb=mp.log(mpq(b))/2;ell=tb-ta
            left,right,off=vf.uppercell_chord(n,hs['lower'],m,a,b)
            M=mp.matrix(4,4)
            def blk(t):return [[mpq(t[0]),mpq(t[1])],[mpq(t[1]),mpq(t[2])]]
            L,Rr,O=blk(left),blk(right),blk(off)
            for p in range(2):
                for q in range(2):
                    M[p,q]=L[p][q];M[2+p,2+q]=Rr[p][q];M[p,2+q]=O[p][q];M[2+p,q]=O[q][p]
            S=fe_condense(n,HL,k,j,ta,ell,N)
            me=mineig(S-M)
            # lower side
            d,cc=vf.lowercell(n,hs['upper'],m,a,b)
            Dn=mp.mpf(vf.D)
            Ml=mp.matrix([[d[0],d[1],cc[0],cc[1]],[d[1],d[2],cc[1],cc[2]],[cc[0],cc[1],d[0],d[1]],[cc[1],cc[2],d[1],d[2]]])/Dn
            El=lin_energy(n,HU,k,j,ta,ell)
            ml=mineig(Ml-El)
            # scale: ell*b
            print(f'n={n} m={m} cell={i} x={float(a):.4g} ell={float(ell):.5f}  min eig(S_FE(gL) - M_upper)={mp.nstr(me,5)}  [rel to ell*x: {mp.nstr(me/(ell*mpq(b)),4)}]   min eig(M_lower - E_lin(gU))={mp.nstr(ml,5)}',flush=True)
            if me <= 0 or ml <= 0:
                raise ArithmeticError(f'numerical cell inequality failed: n={n}, m={m}, cell={i}')
if __name__=='__main__':
    n=int(sys.argv[1]);ms=[int(v) for v in sys.argv[2].split(',')];cells=[int(v) for v in sys.argv[3].split(',')]
    run(n,ms,cells,N=int(sys.argv[4]) if len(sys.argv)>4 else 24)
