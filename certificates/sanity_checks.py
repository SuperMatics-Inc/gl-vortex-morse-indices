"""Non-vacuity checks for the certificates (exit 0 iff every expectation holds).

 (S1) the exact tail numerator N_m / (x^2 (x - x_lo)^2 D_L E_m) agrees with an independent floating-point
      evaluation of V_m - h_tt/h (V_m from its definition, h_tt by finite differences in t) at sample points;
 (S2) a tail parameter set that is numerically bad (m = 3, x_lo = 14.6, s = 2.8, beta = 34) is REJECTED;
 (S3) an upper barrier with too large a constant (D_U(0) = 70000) is REJECTED (b_U changes sign);
 (S4) the exact Rayleigh quotients of the certified lower test vectors agree with the floating-point Ritz
      values recorded by design_lower.py (relative difference < 1e-6);
 (S5) inner supersolutions are not vacuous: extending the m = 2 component past the Bessel zero j_{2,1}^2
      is REJECTED, and the second m = 4 component without the dilation (lam = 1) is REJECTED;
 (S6) the positivity test rejects polynomials with a double root or an interior negative dip.
Usage: python sanity_checks.py [n]   (default n = 4).
"""
import json
import sys
import mpmath as mp
import sympy as sp
from cert_lib import x, positive_on, double_check, Ledger
from vortex import load, potential, tail_poly, bessel_trunc, L_op, barrier_poly
import cert_lower

R = sp.Rational


def rejected(p, a, b=sp.oo):
    try:
        positive_on(p, a, b)
    except ValueError as e:
        return True, str(e)
    return not double_check(p, a, b), 'accepted by Bernstein'


def V_float(n, m, D, xv):
    gsq = xv ** n / D(xv)
    P, W = xv * (1 - gsq), xv * gsq
    A = (n + m) ** 2 - P
    return (n - m) ** 2 - P + A * W / (A + W)


def main(n=4):
    case, DU, DL = load(n)
    led = Ledger()
    Df = sp.lambdify(x, DL)
    for m, xs in ((3, [15.0, 38.3283, 200.0]), (4, [31.0, 80.0]), (2, [55.5, 1000.0])):
        xc, xlo, s, beta = (R(v) for v in case['TAIL'][m])
        Nm = tail_poly(n, m, DL, xlo, s, beta)
        AD, E, Vnum = potential(n, m, DL)
        for xv in xs:
            X = R(str(xv))
            exact = float(Nm.eval(X) / (X ** 2 * (X - xlo) ** 2 * DL.subs(x, X) * E.subs(x, X)))
            mp.mp.dps = 40
            xlo_, s_, beta_ = (mp.mpf(v.p) / v.q for v in (xlo, s, beta))
            logh = lambda t: mp.log(mp.exp(2 * t) - xlo_) - 2 * s_ * t - beta_ * mp.exp(-2 * t)  # noqa: E731
            t0 = mp.log(mp.mpf(xv)) / 2
            d1, d2 = mp.diff(logh, t0, 1), mp.diff(logh, t0, 2)
            fl = V_float(n, m, Df, xv) - float(d2 + d1 ** 2)
            led.report(f'(S1) m={m} x={xv}: exact {exact:.8f} vs float {fl:.8f}',
                       abs(exact - fl) < 1e-8 * max(1, abs(fl)))
    bad, why = rejected(tail_poly(n, 3, DL, R(146, 10), R(28, 10), R(34)), R(367, 25))
    led.report('(S2) bad tail (m=3, x_lo=14.6, s=2.8, beta=34) rejected', bad, f'[{why}]')
    DUbad = sp.expand(DU - DU.subs(x, 0) + 70000)
    bad, why = rejected(-barrier_poly(n, DUbad).as_expr(), 0)
    led.report('(S3) upper barrier with D_U(0)=70000 rejected', bad, f'[{why}]')

    code, rq = cert_lower.run(n)
    data = json.load(open(case['LOWER']))
    led.report('(S4) cert_lower.py passes', code == 0)
    for m, vals in rq.items():
        for i, v in enumerate(vals):
            ev = data[str(m)]['ev'][i]
            led.report(f'(S4) m={m} vector {i}: exact Rayleigh quotient {float(v):.10f} vs Ritz {ev:.10f}',
                       abs(float(v) - ev) < 1e-6 * abs(ev))

    q = bessel_trunc(2, 25)
    bad, why = rejected(q, 0, R(2638, 100))
    led.report('(S5) m=2 inner component extended to x = 26.38 > j_{2,1}^2 rejected', bad, f'[{why}]')
    a, b, sigma, N, lam = case['INNER'][4][1]
    q = bessel_trunc(0, N, 1)
    Lq = sp.Poly(L_op(q, 0), x)
    v = min(mon[0] for mon in Lq.monoms())
    bad_q, why_q = rejected(sigma * q, R(a), R(b))
    bad_L, why_L = rejected(sigma * sp.expand(Lq.as_expr() / x ** v), R(a), R(b))
    led.report('(S5) m=4 second component with lam = 1 rejected', bad_q or bad_L,
               f'[q: {why_q}; L[q]: {why_L}]')

    bad1, why1 = rejected((x - 3) ** 2, 0, 10)
    bad2, why2 = rejected((x - 3) ** 2 * (x - R(31, 10)) ** 2 - R(1, 10 ** 9), 0, 10)
    led.report('(S6) (x-3)^2 on [0,10] rejected', bad1, f'[{why1}]')
    led.report('(S6) (x-3)^2 (x-3.1)^2 - 1e-9 on [0,10] rejected', bad2, f'[{why2}]')
    return led.finish()


if __name__ == '__main__':
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 4))
