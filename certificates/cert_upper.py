"""CERTIFIED (exact rational arithmetic): nu(Q_m) <= K_m for the degree-n blocks, m = 2..2n-2.

Route (del Pino-Liu-Musso-Wei-Yang, arXiv:2608.04976v1, Lemma (slow), Prop. (p:trace), general n):
  Q_m >= Q_{m,g_L} (g_L < f_n, cert_barriers.py) and Q_{m,g_L}[c,d] >= ell_m[d] = int (d_t^2 + V_m d^2) dt,
  t = log r, x = r^2, V_m as in vortex.py with D = D_L.  If R minus K_m points is covered by components
  carrying a positive h with (-d_t^2 + V_m) h > 0, then nu(Q_m) <= K_m.

Checks:
 (U0) D_L > 0, D_L A_m > 0 and E_m = D_L (A_m + W) > 0 on [0, oo); exact rational-function identities
      V_m = Vnum/(D_L E_m) and V_m - (k^2 - x) = Phi_m := W + A_m W/(A_m + W)  (so Phi_m > 0 for x > 0);
      the identity (-d_t^2 + V_m)(x^{k/2} q)/(x^{k/2} q) = Phi_m + L[q]/q for an arbitrary smooth q.
 (U1) inner components (a, b): h = sigma x^{k/2} q, q(x) = p_N^{(k)}(lam x) (truncated Bessel series):
      sigma q > 0 on [a, b] and sigma L[q] / x^v > 0 on [a, b] (x^v the largest power of x dividing L[q]).
      Then (-d_t^2 + V_m) h / h = Phi_m + L[q]/q > 0 on the open component.
 (U2) outer component (x_c, oo): h = (x - x_lo) x^{-s} e^{-beta/x}, x_lo < x_c, and the polynomial
      N_m = (V_m - h_tt/h) x^2 (x - x_lo)^2 D_L E_m is positive on [x_c, oo).
 (U3) the components tile (0, oo) with exactly K_m interior endpoints (the trace points).
Usage: python cert_upper.py [n]   (default n = 4).  Exit code 0 iff every check passes.
"""
import sys
import sympy as sp
from cert_lib import x, Ledger
from vortex import load, potential, potential_symbolic, bessel_trunc, L_op, tail_poly

R = sp.Rational


def check_generic_identity(led):
    """(-d_t^2 + V)(x^{k/2} q)/(x^{k/2} q) - V = -k^2 - (4(k+1) x q' + 4 x^2 q'')/q, with d_t = 2x d_x."""
    k = sp.symbols('k', positive=True)
    X = sp.symbols('X', positive=True)
    q = sp.Function('q')
    h = X ** (k / 2) * q(X)
    Dt = lambda e: 2 * X * sp.diff(e, X)  # noqa: E731
    lhs = sp.simplify(Dt(Dt(h)) / h)
    rhs = k ** 2 + (4 * (k + 1) * X * sp.diff(q(X), X) + 4 * X ** 2 * sp.diff(q(X), X, 2)) / q(X)
    led.report('(U0) h_tt/h = k^2 + (4(k+1) x q\' + 4 x^2 q\'\')/q for h = x^{k/2} q (symbolic, any q)',
               sp.simplify(lhs - rhs) == 0)


def main(n=4):
    case, DU, DL = load(n)
    led = Ledger()
    check_generic_identity(led)
    led.positive('(U0) D_L > 0 on [0,oo)', DL, 0)
    for m in range(2, 2 * n - 1):
        k = abs(n - m)
        AD, E, Vnum = potential(n, m, DL)
        P, W, A, V = potential_symbolic(n, m, DL)
        led.report(f'(U0) m={m}: V_m == Vnum/(D_L E_m)', sp.cancel(V - Vnum / (DL * E)) == 0)
        led.report(f'(U0) m={m}: V_m - (k^2 - x) == W + A_m W/(A_m + W)',
                   sp.cancel(V - (k * k - x) - (W + A * W / (A + W))) == 0)
        led.positive(f'(U0) m={m}: D_L A_m > 0 on [0,oo)', AD, 0)
        led.positive(f'(U0) m={m}: E_m > 0 on [0,oo)', E, 0)

    for m in range(2, 2 * n - 1):
        k = abs(n - m)
        comps = case['INNER'][m]
        xc, xlo, s, beta = (R(v) for v in case['TAIL'][m])
        # (U3) tiling
        ends = [R(comps[0][0])]
        tiled = ends[0] == 0
        for a, b, *_ in comps:
            tiled &= R(a) == ends[-1] and R(b) > R(a)
            ends.append(R(b))
        tiled &= ends[-1] == xc and xlo < xc
        led.report(f'(U3) m={m}: components {[str(e) for e in ends]} + ({xc}, oo); '
                   f'{len(ends) - 1} trace point(s), K_m = {case["K"][m]}',
                   tiled and len(ends) - 1 == case['K'][m])
        # (U1) inner supersolutions
        for a, b, sigma, N, lam in comps:
            a, b, lam = R(a), R(b), R(lam)
            q = bessel_trunc(k, N, lam)
            Lq = sp.Poly(L_op(q, k), x)
            v = min(mon[0] for mon in Lq.monoms())
            Lred = sp.expand(Lq.as_expr() / x ** v)
            tag = f'm={m}, k={k}, sigma={sigma:+d}, N={N}, lam={lam}'
            led.positive(f'(U1) {tag}: sigma q > 0 on [{a}, {b}]', sigma * q, a, b)
            led.positive(f'(U1) {tag}: sigma L[q]/x^{v} > 0 on [{a}, {b}]', sigma * Lred, a, b)
        # (U2) tail
        Nm = tail_poly(n, m, DL, xlo, s, beta)
        led.positive(f'(U2) m={m}: tail N_m > 0 on [{xc}, oo) (x_lo={xlo}, s={s}, beta={beta})',
                     Nm.as_expr(), xc)
    if led.ok:
        for m in range(2, 2 * n - 1):
            print(f'=> nu(Q_{m}) <= {case["K"][m]}')
    return led.finish()


if __name__ == '__main__':
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 4))
