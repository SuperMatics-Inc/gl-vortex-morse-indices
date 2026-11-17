"""Exact positivity tests for polynomials with rational coefficients.

positive_on(p, a, b) certifies p > 0 on the closed interval [a, b] (b may be sympy.oo).
Method: map [a, b] to [0, 1] by x = a + (b - a) z, or [a, oo) by x = a + z/(1 - z) followed by
multiplication with the positive factor (1 - z)^deg; convert to the Bernstein basis of [0, 1]
and subdivide at midpoints (de Casteljau) until every Bernstein coefficient is positive.
Since the Bernstein basis is non-negative and sums to one, the polynomial is then bounded
below on each piece by its smallest coefficient.  All arithmetic is in fractions.Fraction;
no floating-point number enters.  A nonpositive endpoint value is reported as a failure.

double_check(p, a, b) is an independent exact test: sympy's real-root count (Sturm sequence
over the rationals) on [a, b] is zero and p is positive at one rational point of [a, b].
"""
from fractions import Fraction
from math import comb
import sympy as sp

x, z = sp.symbols('x z')


def _to_bernstein(c):
    d = len(c) - 1
    return [sum(Fraction(comb(k, i), comb(d, i)) * c[i] for i in range(k + 1)) for k in range(d + 1)]


def _split(bc):
    """de Casteljau subdivision at 1/2: Bernstein coefficients of the two halves."""
    left, right = [], []
    cur = list(bc)
    while cur:
        left.append(cur[0])
        right.append(cur[-1])
        cur = [(cur[i] + cur[i + 1]) / 2 for i in range(len(cur) - 1)]
    return left, right[::-1]


def _bern_positive(bc, depth=0, maxdepth=60):
    if min(bc) > 0:
        return 1
    if bc[0] <= 0 or bc[-1] <= 0:
        raise ValueError(f'nonpositive endpoint value at depth {depth}')
    if depth >= maxdepth:
        raise ValueError('max subdivision depth reached')
    left, right = _split(bc)
    return _bern_positive(left, depth + 1, maxdepth) + _bern_positive(right, depth + 1, maxdepth)


def _fraction(c):
    c = sp.Rational(c)
    return Fraction(int(c.p), int(c.q))


def positive_on(p, a, b=sp.oo):
    """Return the number of Bernstein pieces used to certify p(x) > 0 on [a, b].

    Raises ValueError when the certificate cannot be produced (in particular when p <= 0
    somewhere on [a, b])."""
    p = sp.Poly(p, x, domain='QQ')
    d = p.degree()
    a = sp.Rational(a)
    if d <= 0:
        if p.LC() > 0:
            return 1
        raise ValueError('nonpositive constant')
    if b is sp.oo:
        expr = sp.cancel(p.as_expr().subs(x, a + z / (1 - z)) * (1 - z) ** d)
    else:
        b = sp.Rational(b)
        assert b > a
        expr = p.as_expr().subs(x, a + z * (b - a))
    pz = sp.Poly(sp.expand(expr), z, domain='QQ')
    c = [_fraction(cf) for cf in reversed(pz.all_coeffs())]
    c += [Fraction(0)] * (d + 1 - len(c))
    return _bern_positive(_to_bernstein(c))


def double_check(p, a, b=sp.oo):
    """Independent exact check: no real root of p in [a, b] and p > 0 at a rational point."""
    p = sp.Poly(p, x, domain='QQ')
    a = sp.Rational(a)
    if p.degree() <= 0:
        return p.LC() > 0
    hi = None if b is sp.oo else sp.Rational(b)
    nroots = p.count_roots(a, hi)
    mid = a + 1 if hi is None else (a + hi) / 2
    return nroots == 0 and p.eval(mid) > 0


class Ledger:
    """Collects PASS/FAIL lines; the calling script exits 0 iff every check passed."""

    def __init__(self):
        self.ok = True
        self.count = 0

    def report(self, name, good, info=''):
        self.ok &= bool(good)
        self.count += 1
        print(f'[{"PASS" if good else "FAIL"}] {name} {info}'.rstrip(), flush=True)
        return bool(good)

    def positive(self, name, p, a, b=sp.oo):
        """Certify p > 0 on [a, b] by Bernstein subdivision and, independently, by Sturm."""
        try:
            boxes = positive_on(p, a, b)
            dc = double_check(p, a, b)
        except ValueError as e:
            boxes, dc = f'FAIL({e})', False
        deg = sp.Poly(p, x).degree()
        return self.report(name, isinstance(boxes, int) and dc,
                           f'(degree {deg}, Bernstein pieces {boxes}, Sturm {dc})')

    def finish(self):
        print(f'{self.count} checks; ' + ('ALL PASS' if self.ok else 'SOME CHECK FAILED'), flush=True)
        return 0 if self.ok else 1
