"""Shared exact (sympy, rational) objects for the degree-n vortex certificates.

Notation (del Pino-Liu-Musso-Wei-Yang, arXiv:2608.04976v1, Sect. 3): x = r^2 = e^{2t},
g^2 = x^n / D(x) for a barrier denominator D,
  P = x (1 - g^2),  W = x g^2,  A_m = (n+m)^2 - P,  V_m = (n-m)^2 - P + A_m W / (A_m + W).
"""
from fractions import Fraction
import sympy as sp
from cases import CASES, asymptotic_coefficients
from cert_lib import x

R = sp.Rational


def load(n):
    c = CASES[n]
    DU = sum(R(v) * x ** i for i, v in enumerate(c['DU']))
    DL = sum(R(v) * x ** i for i, v in enumerate(c['DL']))
    return c, sp.expand(DU), sp.expand(DL)


def barrier_poly(n, D):
    """b(x) with D^2 E_n(g)/g = b(x), g(r) = r^n D(r^2)^{-1/2};
    E_n(g) = g'' + g'/r - n^2 g/r^2 + (1 - g^2) g.  For x > 0, sign E_n(g) = sign b."""
    L = R(n, 2) / x - sp.diff(D, x) / (2 * D)          # d log g / dx
    EoG = 4 * x * (sp.diff(L, x) + L ** 2) + 4 * L - n * n / x + 1 - x ** n / D
    num, den = sp.fraction(sp.together(sp.cancel(EoG * D ** 2)))
    assert sp.Poly(den, x).degree() == 0 and den > 0, den
    return sp.Poly(sp.expand(num / den), x)


def potential(n, m, D):
    """(D*A_m, E_m, Vnum) with E_m = D (A_m + W) and V_m = Vnum / (D E_m)."""
    AD = sp.expand((n + m) ** 2 * D - x * (D - x ** n))
    E = sp.expand((n + m) ** 2 * D - x * D + 2 * x ** (n + 1))
    Vnum = sp.expand((n - m) ** 2 * D * E - x * (D - x ** n) * E + AD * x ** (n + 1))
    return AD, E, Vnum


def potential_symbolic(n, m, D):
    P = x * (D - x ** n) / D
    W = x ** (n + 1) / D
    A = (n + m) ** 2 - P
    return P, W, A, (n - m) ** 2 - P + A * W / (A + W)


def bessel_trunc(k, N, lam=1):
    """p_N^{(k)}(lam x), p_N^{(k)}(y) = sum_{j<=N} (-y/4)^j / (j! (j+k)!).
    J_k(r) = (r/2)^k p_oo^{(k)}(r^2)."""
    return sp.expand(sum(R((-1) ** j, 4 ** j * sp.factorial(j) * sp.factorial(j + k)) * (lam * x) ** j
                         for j in range(N + 1)))


def L_op(q, k):
    """L[q] = -x q - 4(k+1) x q' - 4 x^2 q''; for h = x^{k/2} q one has
    (-d_t^2 + V_m) h / h = Phi_m + L[q]/q with Phi_m = W + A_m W/(A_m + W)."""
    return sp.expand(-x * q - 4 * (k + 1) * x * sp.diff(q, x) - 4 * x ** 2 * sp.diff(q, x, 2))


def tail_poly(n, m, D, xlo, s, beta):
    """N with (V_m - h_tt/h) * x^2 (x - x_lo)^2 D E_m = N(x), h = (x - x_lo) x^{-s} e^{-beta/x}.
    In t: h_t/h = q = 2x/(x - x_lo) - 2s + 2 beta/x and h_tt/h = 2x q_x + q^2."""
    AD, E, Vnum = potential(n, m, D)
    q = 2 * x / (x - xlo) - 2 * s + 2 * beta / x
    hpp = 2 * x * sp.diff(q, x) + q ** 2
    mult = x ** 2 * (x - xlo) ** 2
    num, den = sp.fraction(sp.together(sp.cancel(Vnum * mult - hpp * mult * D * E)))
    assert sp.Poly(den, x).degree() == 0 and den > 0
    return sp.Poly(sp.expand(num / den), x)


def asymptotic_Df(n):
    """Polynomial part (degrees n..n-3) of x^n / f_n^2 from the cited expansion of f_n."""
    a1, a2, a3 = (R(c.numerator, c.denominator) for c in asymptotic_coefficients(n))
    y = sp.symbols('y')
    F = 1 - a1 * y - a2 * y ** 2 - a3 * y ** 3
    inv = sp.series(1 / F ** 2, y, 0, 4).removeO()
    return sp.expand(x ** n * inv.subs(y, 1 / x))


def to_fraction(c):
    c = sp.Rational(c)
    return Fraction(int(c.p), int(c.q))
