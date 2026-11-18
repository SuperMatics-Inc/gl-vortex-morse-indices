"""Check global barriers g_L < f_n < g_U for the degree-n profile with exact rational arithmetic.

    g_U(r)^2 = x^n / D_U(x),  g_L(r)^2 = x^n / D_L(x),  x = r^2   (data in cases.py)

Checks (each positivity by exact Bernstein subdivision and an independent Sturm root count):
 (B0) D_U > 0 and D_L > 0 on [0, oo)  (so g_U, g_L are smooth, g ~ D(0)^{-1/2} r^n at 0, g -> 1).
 (B1) E_n(g_U) < 0 on (0, oo): -b_U > 0, where b = D^2 E_n(g)/g is a polynomial in x.
 (B2) E_n(g_L) > 0 on (0, oo):  b_L > 0.
 (B3) ordering at infinity from the cited expansion f_n = 1 - a1 r^-2 - a2 r^-4 - a3 r^-6 + O(r^-8):
      x^n/f_n^2 = T(x) + O(x^{n-4}), T of degree n with lowest power x^{n-3}.  The leading term of
      T - D_U must be positive (then f_n < g_U for large r) and that of T - D_L negative.
 (B4) D_L A_2 = (n+2)^2 D_L - x (D_L - x^n) > 0 on [0, oo)  (A_m > 0 for every m >= 2).
With the comparison lemma of arXiv:2608.04976v1 (Lemma 2.2) these give g_L < f_n < g_U on (0, oo).
Usage: python cert_barriers.py [n]   (default n = 4).  Exit code 0 iff every check passes.
"""
import sys
import sympy as sp
from cert_lib import x, Ledger
from vortex import load, barrier_poly, asymptotic_Df


def main(n=4):
    case, DU, DL = load(n)
    led = Ledger()
    print(f'n = {n}\nD_U = {DU}\nD_L = {DL}')
    led.positive('(B0) D_U > 0 on [0,oo)', DU, 0)
    led.positive('(B0) D_L > 0 on [0,oo)', DL, 0)
    bU, bL = barrier_poly(n, DU), barrier_poly(n, DL)
    print('b_U(x) =', bU.as_expr())
    print('b_L(x) =', bL.as_expr())
    led.positive('(B1) -b_U > 0 on [0,oo)', -bU.as_expr(), 0)
    led.positive('(B2)  b_L > 0 on [0,oo)', bL.as_expr(), 0)

    T = asymptotic_Df(n)
    print(f'x^{n}/f^2 = {T} + O(x^{n - 4})')
    for name, D, want in (('D_U', DU, 1), ('D_L', DL, -1)):
        top = sp.Poly(sp.expand(T - sum(D.coeff(x, i) * x ** i for i in range(n - 3, n + 1))), x)
        ok = (not top.is_zero) and top.degree() >= n - 3 and sp.sign(top.LC()) == want
        led.report(f'(B3) x^{n}/f^2 - {name} = {top.as_expr()} + O(x^{n - 4})', ok,
                   '(f < g_U near infinity)' if want > 0 else '(f > g_L near infinity)')
    led.positive(f'(B4) (n+2)^2 D_L - x(D_L - x^{n}) > 0 on [0,oo)',
                 sp.expand((n + 2) ** 2 * DL - x * (DL - x ** n)), 0)

    lo, hi = 1 / float(DL.subs(x, 0)) ** 0.5, 1 / float(DU.subs(x, 0)) ** 0.5
    print(f'consequence: kappa_{n} = lim f/r^{n} lies in [D_L(0)^-1/2, D_U(0)^-1/2] = [{lo:.6f}, {hi:.6f}]')
    return led.finish()


if __name__ == '__main__':
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 4))
