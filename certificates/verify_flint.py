"""Check the degree-four certificates with FLINT rationals and a separate Sturm implementation.
No imports from the original certificates, SymPy, NumPy or SciPy. All polynomial
signs include both finite endpoints; infinity is handled by leading coefficients.
JSON binary64 vectors only specify exact rational inputs, as in the manuscript.
"""
import json
from pathlib import Path
from fractions import Fraction
from math import factorial
from flint import fmpq as Q, fmpq_poly as Poly
import flint

x = Poly([0, 1])
ONE = Poly([1])
DU = Poly([55000, 6128, 288, 16, 1])
DL = Poly([500000, 6160, 288, 16, 1])
COUNTS = {2: 1, 3: 1, 4: 2, 5: 1, 6: 1}
INNER = {2: [(0,'2637/100',1,25,'1')], 3: [(0,'367/25',1,21,'1')],
         4: [(0,'144579/25000',1,15,'1'), ('144579/25000','3047/100',-1,20,'100001/100000')],
         5: [(0,'367/25',1,21,'1')], 6: [(0,'2637/100',1,25,'1')]}
TAIL = {2: ('2637/100','26','7/4','26'), 3: ('367/25','1467/100','281/100','339/10'),
        4: ('3047/100','152/5','14/5','70'), 5: ('367/25','73/5','13/4','40'),
        6: ('2637/100','26','2','30')}
checks = 0

def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(label)
    print('PASS', label, flush=True)

def variations(signs):
    nz = [s for s in signs if s]
    return sum(a != b for a, b in zip(nz, nz[1:]))

def sign(v):
    return (v > 0) - (v < 0)

def positive(p, lo=0, hi=None):
    """Exact Sturm count, including repeated roots and strict endpoint checks."""
    lo = Q(lo)
    if not p or p(lo) <= 0 or (hi is not None and p(Q(hi)) <= 0):
        return False
    seq = [p]
    if p.degree() > 0:
        seq.append(p.derivative())
        while seq[-1]:
            rem = -(seq[-2] % seq[-1])
            if not rem:
                break
            # Positive scaling contains coefficient growth, preserves Sturm signs.
            rem /= abs(rem.leading_coefficient())
            seq.append(rem)
    left = variations([sign(t(lo)) for t in seq])
    right = variations([sign(t.leading_coefficient() if hi is None else t(Q(hi))) for t in seq])
    return left == right

def residual(D):
    return D*D-x**4*D-2*x*D*D.derivative().derivative()+3*x*D.derivative()**2-10*D*D.derivative()

def potential(m):
    ad = (4+m)**2*DL-x*(DL-x**4)
    e = ad+x**5
    vn = (4-m)**2*DL*e-x*(DL-x**4)*e+ad*x**5
    return ad, e, vn

def trunc(k, n, lam):
    return Poly([(-Q(lam)/4)**j / (factorial(j)*factorial(j+k)) for j in range(n+1)])

def tail(m, low, s, beta):
    ad, e, vn = potential(m)
    low, s, beta = Q(low), Q(s), Q(beta)
    den = x*(x-low)
    num = 2*x*x-2*s*den+2*beta*(x-low)
    # q=num/den; 2x q' + q^2, clearing den^2.
    hpp_num = 2*x*(num.derivative()*den-num*den.derivative())+num*num
    return vn*den*den-hpp_num*DL*e

def integral(p, a, b):
    prim = p.integral()
    return prim(Q(b))-prim(Q(a))

def tests(entry, count):
    X, degree = Q(entry['X']), entry['p']
    z = 2*x/X-1
    legs = [ONE, z]
    for i in range(1, degree):
        legs.append(((2*i+1)*z*legs[-1]-i*legs[-2])/(i+1))
    out = []
    for col in range(count):
        v = [row[col] for row in entry['vec']]
        big = max(abs(c) for c in v)
        fractions = [Fraction(c/big).limit_denominator(1000000) for c in v]
        c = [Q(t.numerator,t.denominator) for t in fractions]
        cut = (1-x/X)**2
        A = cut*sum((c[i]*legs[i] for i in range(degree+1)), Poly([]))
        B = cut*sum((c[degree+1+i]*legs[i] for i in range(degree+1)), Poly([]))
        require(A(X)==0 and A.derivative()(X)==0 and B(X)==0 and B.derivative()(X)==0, 'double support cutoff')
        out.append((A, B))
    return out

def form(m, v, w, X, pieces):
    k, j = abs(4-m), abs(4-m)+2
    A, B = v
    C, D = w
    # Derivatives of r^k A(r^2), expanded independently of original script.
    if k:
        kinetic_a = x**(k-1)*(k*A+2*x*A.derivative())*(k*C+2*x*C.derivative())
        cent_a = k*k*x**(k-1)*A*C
    else:
        kinetic_a = 4*x*A.derivative()*C.derivative()
        cent_a = Poly([])
    kinetic_b = x**(j-1)*(j*B+2*x*B.derivative())*(j*D+2*x*D.derivative())
    constant = (kinetic_a+cent_a+kinetic_b+(4+m)**2*x**(j-1)*B*D-x**k*A*C-x**j*B*D)/2
    coupling = x**k*A*C+x**j*B*D+x**(k+1)*(A*D+C*B)/2
    return integral(constant,0,X)+sum(integral(U*coupling,a,b) for a,b,U in pieces)

def main():
    print('Independent python-flint', flint.__version__)
    for label, p in [('DU',DU),('DL',DL),('-bU',-residual(DU)),('bL',residual(DL))]:
        require(positive(p), label+' strictly positive on [0,infinity)')
    # Formal inverse-square series algebra; rigorous remainder remains external.
    f = Poly([1,-8,-48,-896])
    inv = Poly([1,16,288,6144])
    require(all((f*f*inv)[i] == (1 if i==0 else 0) for i in range(4)), 'asymptotic inverse series through order 3')
    require(Q(6144)-DU[1]>0 and Q(6144)-DL[1]<0, 'strict ordering at infinity conditional on cited remainder')
    for m in COUNTS:
        ad,e,vn = potential(m)
        require(positive(ad) and positive(e), f'm={m} positive elimination denominators')
        # Clearing denominators in V-(k^2-x) = W+AW/(A+W).
        require(vn-((4-m)**2-x)*DL*e == x**5*e+ad*x**5, f'm={m} scalar identity')
        last = Q(0)
        for a,b,sigma,N,lam in INNER[m]:
            a,b = Q(a),Q(b)
            require(a==last and b>a, f'm={m} interval adjacency')
            last = b
            q = trunc(abs(4-m),N,lam)
            L = -x*q-4*(abs(4-m)+1)*x*q.derivative()-4*x*x*q.derivative().derivative()
            valuation = next(i for i,c in enumerate(L.coeffs()) if c)
            reduced = Poly(L.coeffs()[valuation:])
            require(positive(sigma*q,a,b), f'm={m} inner polynomial sign N={N}')
            require(positive(sigma*reduced,a,b), f'm={m} inner strict residual N={N}')
        xc,low,s,beta = TAIL[m]
        require(last==Q(xc) and Q(low)<last and len(INNER[m])==COUNTS[m], f'm={m} full component coverage')
        require(positive(tail(m,low,s,beta),xc), f'm={m} tail residual to infinity')
    saved = json.loads(Path(__file__).with_name('lower_tests_n4.json').read_text())
    pieces=[]
    last=Q(0)
    for entry in saved['majorant']:
        a,b=Q(entry['lo']),Q(entry['hi'])
        U=Poly([Q(c) for c in entry['coeffs']])
        require(a==last and b>a and positive(U*DU-x**4,a,b), f'majorant [{a},{b}]')
        pieces.append((a,b,U)); last=b
    for m,count in COUNTS.items():
        entry=saved[str(m)]; X=Q(entry['X'])
        require(X==last, f'm={m} support exactly covered')
        vs=tests(entry,count)
        M=[[form(m,v,w,X,pieces) for w in vs] for v in vs]
        require(M[0][0]<0, f'm={m} exact M11 negative: {float(M[0][0]):.12g}')
        if count==2:
            det=M[0][0]*M[1][1]-M[0][1]*M[1][0]
            require(M[0][1]==M[1][0] and det>0, f'm=4 symmetric negative definite, det={float(det):.12g}')
            print('EXACT m=4 determinant', det)
        print(f'EXACT m={m} M11',M[0][0])
    # Adversarial regression checks; no original checker imported.
    require(not positive(-residual(DU+15000)), 'REJECT upper barrier constant 70000')
    require(not positive(tail(3,'146/10','28/10','34'),'367/25'), 'REJECT bad tail')
    require(not positive(trunc(2,25,1),0,'2638/100'), 'REJECT extended inner interval')
    require(not positive(-trunc(0,20,1),'144579/25000','3047/100'), 'REJECT undilated central interval')
    require(not positive((x-3)**2,0,10), 'REJECT interior double root')
    require(not positive((x-3)**2*(x-Q('31/10'))**2-Q('1/1000000000'),0,10), 'REJECT narrow negative dip')
    require(not positive(-ONE,0,160), 'REJECT negative majorant')
    v=tests(saved['4'],1)[0]
    q=form(4,v,v,Q(160),pieces)
    require(q*q-q*q==0, 'REJECT duplicate two-vector basis: determinant zero')
    require(2*sum(COUNTS.values())==12, 'real index arithmetic = 12')
    import mpmath as mp
    mp.iv.dps=40
    print('INTERVAL diagnostic only: kappa enclosure endpoints',1/mp.iv.sqrt(500000),1/mp.iv.sqrt(55000))
    print(f'{checks} checks; ALL PASS; active signs use exact FLINT rationals, not intervals')

if __name__=='__main__':
    main()
