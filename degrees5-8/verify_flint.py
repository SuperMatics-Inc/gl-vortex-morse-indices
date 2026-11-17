"""Independent verifier: FLINT polynomials and rationals only, no designer imports.
Rejects malformed profile identities, upper inertia witnesses, lower trial/Gram witnesses.
"""
import json,sys,copy,time,gzip
from pathlib import Path
from math import comb
from flint import fmpq as Q,fmpq_poly as P,fmpz_mat
R=10**15;D=10**12
def fl(q,d=R):return Q(int(q.numerator)*d//int(q.denominator),d)
def ce(q,d=R):return -fl(-q,d)
def wholeup(q,d):return -((-int(q.numerator)*d)//int(q.denominator))
def near(q,d):
    # Python round has ties-to-even. Exact integer implementation.
    a=int(q.numerator)*d;b=int(q.denominator);v,r=divmod(a,b)
    return v+(2*r>b or (2*r==b and v%2!=0))
def positive(poly,a=Q(0),b=Q(1)):
    if not poly:raise ValueError('zero polynomial')
    lin=P([a,b-a]);p=P([])
    for c in reversed(poly.coeffs()):p=p*lin+c
    deg=p.degree();co=p.coeffs();co+= [Q(0)]*(deg+1-len(co))
    bs=[sum((co[j]*Q(comb(k,j),comb(deg,j)) for j in range(k+1)),Q(0)) for k in range(deg+1)]
    stack=[(bs,0)];pieces=0
    while stack:
        v,d=stack.pop()
        if min(v)>0:pieces+=1;continue
        if v[0]<=0 or v[-1]<=0 or d>=24:raise ValueError('nonpositive/unresolved polynomial')
        levels=[v]
        while len(levels[-1])>1:levels.append([(u+w)/2 for u,w in zip(levels[-1][:-1],levels[-1][1:])])
        stack.extend([([row[0] for row in levels],d+1),([row[-1] for row in reversed(levels)],d+1)])
    return pieces
EXPECTED_INDICES={5:20,6:30,7:36,8:48}

def checked_barrier(n,H,side):
    """Check the mathematical side, not a caller-supplied direction label."""
    assert side in ('lower','upper') and type(n) is int and n>=1
    y=P([0,1]);sgn=1 if side=='upper' else -1
    assert H(Q(1))==1 and H.derivative()(Q(1))==Q(-(n-1),2)
    B=-n*(2*(n+1)-(n+2)*y)*H+4*(1-y)*((n+1)-(n+2)*y)*H.derivative()+4*y*(1-y)**2*H.derivative().derivative()
    red,rem=divmod(-sgn*(n*n*H*(1-y**n*H*H)+(1-y)*B),(1-y)**2)
    assert not rem
    positive(H);positive(red);positive(n*H+2*y*H.derivative())
    cg=n**4*(H.derivative().derivative()(Q(1))+Q(5-n*n,4))
    assert n*n*red(Q(1))==sgn*(cg+2*n*n)>0
    return red,cg

def matching_counts(n,lo,hi,lc,uc):
    required=set(map(str,range(2,2*n-1)))
    assert n in EXPECTED_INDICES and lo==hi==EXPECTED_INDICES[n]
    assert set(lc)==set(uc)==required and lc==uc


def profiles(rec):
    n=rec['n'];y=P([0,1]);hs={}
    for name,sgn in [('upper',1),('lower',-1)]:
        row=rec['barriers'][name];H=P([Q(v) for v in row['H']]);hs[name]=H
        red,cg=checked_barrier(n,H,name)
        assert red==P([Q(v) for v in row['residual_reduced']])
    delta,rem=divmod(hs['upper']-hs['lower'],(1-y)**2);assert not rem;positive(delta)
    return hs

def G(n,H,x):
    y=x/(n*n+x);return y**n*H(y)**2

def length(a,b):
    z=(b-a)/(b+a);s=z+z**3/3+z**5/5
    return fl(s),ce(s+z**7/(7*(1-z*z)))

def lowercell(n,H,m,a,b):
    lo,hi=length(a,b);w=ce(b*G(n,H,b));k=abs(n-m);j=n+m
    aa=k*k-a+2*w;bb=j*j-a+2*w;cc=w
    C=max(abs(aa)+abs(cc),abs(bb)+abs(cc));err=ce((hi-lo)/hi*C)
    aa+=err;bb+=err
    return (near(1/lo+hi*aa/3,D)+4,near(hi*cc/3,D),near(1/lo+hi*bb/3,D)+4),(near(-1/lo+hi*aa/6,D),near(hi*cc/6,D),near(-1/lo+hi*bb/6,D))

def lower(rec,H):
    n=rec['degree'];xs=list(map(Q,rec['x_grid']));assert xs[0]>0 and all(a<b for a,b in zip(xs[:-1],xs[1:]));assert rec['matrix_denominator']==D and rec['vector_denominator']==10**7
    checked_barrier(n,H,'upper')
    assert set(rec['channels'])==set(map(str,range(2,2*n-1)))
    count=0
    for sm,r in rec['channels'].items():
        m=int(sm);K=r['count'];assert 1<=K<=3
        V=r['vectors'];assert len(V)==2*len(xs) and all(len(v)==K for v in V) and V[-1]==V[-2]==[0]*K
        assert all(type(t) is int for row in V for t in row)
        A=[[0]*K for _ in range(K)]
        k=abs(n-m);j=n+m;x=xs[0];g=G(n,H,x)
        org=[wholeup(k+3*x*g/(2*(k+1)),D),wholeup(j+3*x*g/(2*(j+1)),D)]
        for p in range(K):
            for q in range(K):A[p][q]+=org[0]*V[0][p]*V[0][q]+org[1]*V[1][p]*V[1][q]
        for i,(a,b) in enumerate(zip(xs[:-1],xs[1:])):
            d,c=lowercell(n,H,m,a,b)
            M=[[d[0],d[1],c[0],c[1]],[d[1],d[2],c[1],c[2]],[c[0],c[1],d[0],d[1]],[c[1],c[2],d[1],d[2]]]
            v=V[2*i:2*i+4]
            for p in range(K):
                for q in range(K):A[p][q]+=sum(M[s][t]*v[s][p]*v[t][q] for s in range(4) for t in range(4))
        assert A==[[int(v) for v in row] for row in r['gram_integer']]
        for sz in range(1,K+1):
            det=int(fmpz_mat([row[:sz] for row in A[:sz]]).det());assert det==int(r['leading_minors'][sz-1]) and (-1)**sz*det>0
        count+=K
    assert rec['lower_index']==2*count
    return 2*count

# Interval routines implemented as endpoint tuples, distinct from the class-based designer.
def plus(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def times(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def scale(a,s):return times(a,(s,s))
def sub(a,b):return plus(a,neg(b))
def inv(a):
    assert not a[0]<=0<=a[1];return (1/a[1],1/a[0])
def uppercell(n,H,m,a,b):
    lo,h=length(a,b);w=fl(a*G(n,H,a));k=abs(n-m);j=n+m
    aa=k*k-b+2*w;bb=j*j-b+2*w;cc=w;C=max(abs(aa)+abs(cc),abs(bb)+abs(cc));assert 9-C*h*h>0
    eta=ce(C*C*h*h/(9-C*h*h));err=ce((h-lo)/h*(C+eta));aa-=eta+err;bb-=eta+err
    return (1/h+h*aa/3,h*cc/3,1/h+h*bb/3),(-1/h+h*aa/6,h*cc/6,-1/h+h*bb/6)
def uppercell_chord(n,H,m,a,b):
    lo,h=length(a,b);k=abs(n-m);j=n+m
    ws=[fl(x*G(n,H,x)) for x in [a,b]]
    ends=[(k*k-x+2*w,w,j*j-x+2*w) for x,w in zip([a,b],ws)]
    A,B=ends;delta=ce(h*h*b*31/2)
    C=max(abs(v[i])+abs(v[1]) for v in ends for i in [0,2])+delta
    assert 9-C*h*h>0
    eta=ce(C*C*h*h/(9-C*h*h));loss=delta+eta+ce((h-lo)/h*(C+eta))
    left=[];right=[];off=[]
    for i in range(3):
        diag=i in [0,2]
        left.append((1/h if diag else 0)+h*(A[i]/4+B[i]/12-(loss/3 if diag else 0)))
        right.append((1/h if diag else 0)+h*(A[i]/12+B[i]/4-(loss/3 if diag else 0)))
        off.append((-1/h if diag else 0)+h*((A[i]+B[i])/12-(loss/6 if diag else 0)))
    return left,right,off
def upper(rec,H):
    n=rec['degree'];xs=list(map(Q,rec['x_grid']));assert xs[0]>0 and all(a<b for a,b in zip(xs[:-1],xs[1:]));assert xs[0]<=1
    checked_barrier(n,H,'lower')
    den=int(rec['precision_denominator']);count=0;y=P([0,1]);X=xs[-1]
    assert den>0 and int(rec['round_denominator'])==R
    assert rec.get('scheme') in [None,'chord_v1'] and rec.get('inverse_scheme') in [None,'schur_v1']
    chord=rec.get('scheme')=='chord_v1'
    if chord:
        g=y**n*H*H
        z=g+y*(1-y)*(3-2*y)*g.derivative()+y*y*(1-y)**2*g.derivative().derivative()
        positive(10-z);positive(10+z)
    assert set(rec['channels'])==set(map(str,range(2,2*n-1)))
    for sm,r in rec['channels'].items():
        m=int(sm);k=abs(n-m);j=n+m;W=n*n*y**(n+1)*H*H
        robin=r['tail'].get('robin',0);assert type(robin) is int and robin>=0
        aa=(k*k-robin*robin)*(1-y)-n*n*y+2*W;bb=(j*j-robin*robin)*(1-y)-n*n*y+2*W
        dt,rem=divmod(aa*bb-W*W,1-y);assert not rem;positive(aa,X/(n*n+X));positive(dt,X/(n*n+X))
        build=uppercell_chord if chord else uppercell
        cs=[build(n,H,m,a,b) for a,b in zip(xs[:-1],xs[1:])];prev=None;ct=0
        assert len(r['pivots'])==len(xs)
        for i,piv in enumerate(r['pivots']):
            d=[Q(0)]*3
            if i:
                d=[u+v for u,v in zip(d,cs[i-1][1 if chord else 0])]
            if i<len(cs):d=[u+v for u,v in zip(d,cs[i][0])]
            if not i:d[0]+=k-xs[0]/(k+1);d[2]+=j-xs[0]/(j+1)
            if i==len(cs):d[0]+=robin;d[2]+=robin
            expected=[(v,v) for v in d]
            if i:
                a,b,c=prev
                if rec.get('inverse_scheme')=='schur_v1' and not a[0]<=0<=a[1] and max(abs(a[0]),abs(a[1]))>=max(abs(c[0]),abs(c[1])):
                    ia=inv(a);v=times(b,ia);z=inv(sub(c,times(b,v)));q=times(neg(v),z);p=plus(ia,times(times(v,v),z))
                elif rec.get('inverse_scheme')=='schur_v1' and not c[0]<=0<=c[1]:
                    ic=inv(c);v=times(b,ic);p=inv(sub(a,times(b,v)));q=times(neg(v),p);z=plus(ic,times(times(v,v),p))
                else:
                    di=inv(sub(times(a,c),times(b,b)));p,q,z=times(c,di),neg(times(b,di)),times(a,di)
                u,v,w=cs[i-1][2 if chord else 1]
                terms=[plus(plus(scale(p,u*u),scale(q,2*u*v)),scale(z,v*v)),plus(plus(scale(p,u*v),scale(q,u*w+v*v)),scale(z,v*w)),plus(plus(scale(p,v*v),scale(q,2*v*w)),scale(z,w*w))]
                expected=[sub(a,b) for a,b in zip(expected,terms)]
            got=[(Q(int(v[0]),den),Q(int(v[1]),den)) for v in piv]
            assert all(g[0]<=e[0]<=e[1]<=g[1] for g,e in zip(got,expected))
            a,b,c=got;dt=sub(times(a,c),times(b,b))
            if dt[1]<0:ct+=1
            elif dt[0]>0 and a[1]<0:ct+=2
            else:assert dt[0]>0 and a[0]>0
            prev=got
        assert ct==r['count'];count+=ct
    assert rec['upper_index']==2*count;return 2*count

def main():
    root=Path(__file__).parent;cert=Path(sys.argv[1] if len(sys.argv)>1 else root/'certificates');results={};mut=[]
    files=sorted((root/'profiles').glob('profile_n*.json'))
    assert {json.loads(p.read_text())['n'] for p in files}==set(EXPECTED_INDICES) and len(files)==len(EXPECTED_INDICES)
    for p in files:
        r=json.loads(p.read_text());hs=profiles(r);n=r['n'];results[n]={'profile':'PASS'};accepted={}
        bad=copy.deepcopy(r);bad['barriers']['upper']['residual_reduced'][0]='-1'
        try:profiles(bad)
        except (AssertionError,ValueError):mut.append(f'n{n}: corrupted residual rejected')
        else:raise AssertionError('residual mutant accepted')
        for kind,fun in [('lower',lower),('upper',upper)]:
            candidates=sorted(cert.glob(f'{kind}_n{n}*.json*'));assert len(candidates)==1
            for f in candidates:
                c=json.loads(gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes());assert c['degree']==n;H=hs['upper' if kind=='lower' else 'lower'];value=fun(c,H);results[n][f.name]=value
                accepted[kind]=(value,{m:v['count'] for m,v in c['channels'].items()})
                bad=copy.deepcopy(c)
                if kind=='lower':bad['channels']['2']['vectors'][0][0]+=12345
                else:bad['channels']['2']['pivots'][0][0]=['0','0']
                try:fun(bad,H)
                except (AssertionError,ValueError):mut.append(f'{f.name}: corrupted {kind} witness rejected')
                else:raise AssertionError('certificate mutant accepted')
                print('FLINT VERIFIED',f.name,value,flush=True)
                missing=copy.deepcopy(c);del missing['channels']['2']
                try:fun(missing,H)
                except (AssertionError,ValueError):mut.append(f'{f.name}: missing angular channel rejected')
                else:raise AssertionError('missing-channel mutant accepted')
        matching_counts(n,accepted['lower'][0],accepted['upper'][0],accepted['lower'][1],accepted['upper'][1])
        results[n]['matched_expected_index']=EXPECTED_INDICES[n]
    out=Path(__import__('os').environ.get('SM_OUT','.'));(out/'flint_verification.json').write_text(json.dumps({'results':results,'mutants':mut},indent=2));print('FLINT ALL PASS',results,'MUTANTS',len(mut),flush=True)
if __name__=='__main__':main()
