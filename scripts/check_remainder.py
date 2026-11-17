"""Exact symbolic checks for the exterior profile-remainder lemma.

This checks algebra, not the analytic maximum-principle argument.
Run: python -B scripts/check_remainder.py
"""
import json
import sympy as sp

z, s = sp.symbols('z s')
a = -s / 2
b = -s * (s + 8) / 8
d = -s * (s**2 + 32*s + 128) / 16
F = 1 + a*z + b*z**2 + d*z**3
residual = sp.Poly(sp.expand(4*z**3*sp.diff(F, z, 2) + 4*z**2*sp.diff(F, z) - s*z*F + (1-F**2)*F), z)
assert all(residual.nth(j) == 0 for j in range(4))
assert residual.degree() <= 9
assert [F.coeff(z, j).subs(s, 16) for j in range(4)] == [1, -8, -48, -896]
C, r = sp.symbols('C r', positive=True)
h = C*r**-8
assert sp.simplify(-sp.diff(h, r, 2)-sp.diff(h, r)/r + 64*C*r**-10) == 0
assert sp.Rational(64, 128) == sp.Rational(1, 2)
print(json.dumps({'status': 'EXACT ALGEBRA PASS', 'cancelled_powers': [0, 1, 2, 3],
                  'residual_coefficients': {str(j): str(sp.factor(residual.nth(j))) for j in range(4, 10)},
                  'degree_four_coefficients': [1, -8, -48, -896],
                  'comparison_derivative_factor': 64, 'R_squared_threshold': 128}, indent=2))
