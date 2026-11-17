"""CERTIFIED (exact rational arithmetic): ind(Q_m) >= K_m for the degree-n blocks, m = 2..2n-2.

Monotonicity: for G >= f_n^2 on [0, X], Q_{m,G} - Q_m = 2 int (G - f_n^2)(a^2 + ab + b^2) r dr >= 0 on
pairs supported in r^2 <= X (a^2 + ab + b^2 >= 0 pointwise).  So a subspace on which Q_{m,G} is
negative definite is a negative subspace of Q_m.  We take for G the piecewise polynomial U saved in
the test-data file (rational coefficients; how it was found plays no role) and check
  (L0) the pieces tile [0, x_end] and U_i(x) D_U(x) - x^n > 0 on each [x_i, x_{i+1}]  (exact Bernstein
       + Sturm), so U > g_U^2 > f_n^2 there (cert_barriers.py gives f_n < g_U).
Test pairs (saved Ritz vectors, coefficients rationalised): a = r^k A(x), b = r^{k+2} B(x),
x = r^2, k = |n-m|, A, B polynomials with a double zero at x = X <= x_end (extended by 0; in the form
domain).  In x,  Q_{m,U}(v, w) = int_0^X Poly_vw dx + sum_i int_{x_i}^{x_{i+1}} U_i Wg_vw dx  (exact),
  Poly = [x^{k-1}(kA+2xA')^2 + k^2 x^{k-1}A^2 - x^k A^2 + x^{k+1}(jB+2xB')^2 + (n+m)^2 x^{k+1} B^2
          - x^{k+2} B^2]/2  (j = k+2; for k = 0 the a-part is 2x A'^2 - A^2/2),
  Wg   = x^k A^2 + x^{k+2} B^2 + x^{k+1} A B   (bilinear versions by polarisation).
 (L1) K_m = 1: Q_{m,U}[v] < 0.  K_m = 2: M_11 < 0 and det M > 0 (M the Gram matrix of Q_{m,U}).
Usage: python cert_lower.py [n]   (default n = 4).  Exit code 0 iff every check passes.
"""
import json
import sys
from fractions import Fraction
import sympy as sp
from cert_lib import x, Ledger, positive_on, double_check
from vortex import load

R = sp.Rational


def rat(v, den=10 ** 6):
    f = Fraction(v).limit_denominator(den)
    return R(f.numerator, f.denominator)


def build_tests(entry, K):
    """Rationalised test polynomials (A, B) for the K saved Ritz vectors."""
    X, p = entry['X'], entry['p']
    s = x / X
    cut = (1 - s) ** 2
    nb = p + 1
    tests = []
    for col in range(K):
        vec = [row[col] for row in entry['vec']]
        big = max(abs(c) for c in vec)
        co = [rat(c / big) for c in vec]
        A = sp.expand(cut * sum(co[i] * sp.legendre(i, 2 * s - 1) for i in range(nb)))
        B = sp.expand(cut * sum(co[nb + i] * sp.legendre(i, 2 * s - 1) for i in range(nb)))
        tests.append((sp.Poly(A, x, domain='QQ'), sp.Poly(B, x, domain='QQ')))
    return tests


def majorant(n, DU, saved, led):
    """Saved pieces (x_i, x_{i+1}, U_i); certifies tiling from 0 and U_i D_U - x^n > 0 on each piece."""
    pieces, end = [], R(0)
    for piece in saved:
        xa, xb = R(piece['lo']), R(piece['hi'])
        U = sp.Poly(sum(R(c) * x ** i for i, c in enumerate(piece['coeffs'])), x, domain='QQ')
        test = sp.expand(U.as_expr() * DU - x ** n)
        try:
            nb = positive_on(test, xa, xb)
            ok = xa == end and xb > xa and double_check(test, xa, xb)
        except ValueError as err:
            nb, ok = str(err), False
        led.report(f'(L0) U D_U - x^{n} > 0 on [{xa}, {xb}]', ok,
                   f'(degree {U.degree()}, Bernstein pieces {nb}, Sturm {ok})')
        if not ok:
            return None, end
        pieces.append((xa, xb, U))
        end = xb
    return pieces, end


def definite(p, a, b):
    P = p.integrate()
    return P.eval(b) - P.eval(a)


def bilinear(n, m, v, w, X, pieces):
    k = abs(n - m)
    j = k + 2
    J2 = (n + m) ** 2
    (A1, B1), (A2, B2) = v, w
    dA1, dA2, dB1, dB2 = A1.diff(x), A2.diff(x), B1.diff(x), B2.diff(x)
    X1 = sp.Poly(x, x, domain='QQ')
    half = R(1, 2)
    if k == 0:
        aa = 2 * X1 * dA1 * dA2 - A1 * A2 * half
    else:
        aa = (X1 ** (k - 1) * (k * A1 + 2 * X1 * dA1) * (k * A2 + 2 * X1 * dA2)
              + k * k * X1 ** (k - 1) * A1 * A2 - X1 ** k * A1 * A2) * half
    bb = (X1 ** (j - 1) * (j * B1 + 2 * X1 * dB1) * (j * B2 + 2 * X1 * dB2)
          + J2 * X1 ** (j - 1) * B1 * B2 - X1 ** j * B1 * B2) * half
    Wg = X1 ** k * A1 * A2 + X1 ** j * B1 * B2 + X1 ** (k + 1) * (A1 * B2 + A2 * B1) * half
    val = definite(aa + bb, 0, X)
    for xa, xb, U in pieces:
        val += definite(U * Wg, xa, xb)
    return val


def mass(n, m, v, X):
    k = abs(n - m)
    A, B = v
    X1 = sp.Poly(x, x, domain='QQ')
    return definite((X1 ** k * A * A + X1 ** (k + 2) * B * B) * R(1, 2), 0, X)


def run(n=4):
    """Returns (exit code, {m: exact Rayleigh quotients of the test vectors})."""
    case, DU, DL = load(n)
    led = Ledger()
    data = json.load(open(case['LOWER']))
    pieces, end = majorant(n, DU, data['majorant'], led)
    if pieces is None:
        return led.finish(), {}
    results = {}
    for m in range(2, 2 * n - 1):
        K = case['K'][m]
        entry = data[str(m)]
        X = R(entry['X'])
        if not 0 < X <= end:
            led.report(f'(L1) m={m}: the majorant covers [0, X = {X}]', False)
            continue
        use = [(a, min(b, X), U) for a, b, U in pieces if a < X]
        tests = build_tests(entry, K)
        M = [[None] * K for _ in range(K)]
        for i in range(K):
            for jj in range(i, K):
                M[i][jj] = M[jj][i] = bilinear(n, m, tests[i], tests[jj], X, use)
        if K == 1:
            led.report(f'(L1) m={m}: Q_(m,U)[v] = {float(M[0][0]):.10g} < 0', M[0][0] < 0, f'(X={X}, p={entry["p"]})')
        else:
            assert K == 2
            det = M[0][0] * M[1][1] - M[0][1] ** 2
            led.report(f'(L1) m={m}: M11 = {float(M[0][0]):.10g} < 0, M22 = {float(M[1][1]):.10g}, '
                       f'M12 = {float(M[0][1]):.4g}, det = {float(det):.10g} > 0',
                       M[0][0] < 0 and det > 0, f'(X={X}, p={entry["p"]})')
        rq = [M[i][i] / mass(n, m, tests[i], X) for i in range(K)]
        results[m] = rq
        print(f'      m={m}: Rayleigh quotients Q_(m,U)[v_i]/|v_i|^2 = {[f"{float(v):.9f}" for v in rq]}; '
              f'float Ritz values (g_U) {[round(v, 9) for v in entry["ev"][:K]]}')
    if led.ok:
        for m in range(2, 2 * n - 1):
            print(f'=> ind(Q_{m}) >= {case["K"][m]}')
    return led.finish(), results


if __name__ == '__main__':
    sys.exit(run(int(sys.argv[1]) if len(sys.argv) > 1 else 4)[0])
