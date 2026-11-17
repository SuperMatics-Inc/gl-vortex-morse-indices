"""Certificate data for the degree-n vortex (all numbers exact rationals, given as strings).

DU, DL: barrier denominators, low -> high coefficients; g^2 = x^n / D, x = r^2.
K:      claimed index of the radial block Q_m, m = 2 .. 2n-2.
INNER:  for each m, the inner components of the scalar supersolution argument, left to right:
        (a, b, sigma, N, lam) means h = sigma * x^{k/2} p_N^{(k)}(lam * x) on the x-interval (a, b),
        p_N^{(k)}(y) = sum_{j<=N} (-y/4)^j / (j! (j+k)!),  k = |n - m|.
TAIL:   for each m, (x_c, x_lo, s, beta): h = (x - x_lo) x^{-s} exp(-beta/x) on (x_c, oo).
        The trace points are the interior endpoints, so there are K_m of them.
LOWER:  JSON file with the rational majorant U and the negative test subspaces (from design_lower.py).
LOWER_BREAKS, LOWER_DEG: pieces and degree of the polynomial majorant U >= x^n / D_U (design_lower.py input).
"""

CASES = {
    4: dict(
        DU=['55000', '6128', '288', '16', '1'],
        DL=['500000', '6160', '288', '16', '1'],
        K={2: 1, 3: 1, 4: 2, 5: 1, 6: 1},
        INNER={
            2: [('0', '2637/100', 1, 25, '1')],
            3: [('0', '367/25', 1, 21, '1')],
            4: [('0', '144579/25000', 1, 15, '1'),
                ('144579/25000', '3047/100', -1, 20, '100001/100000')],
            5: [('0', '367/25', 1, 21, '1')],
            6: [('0', '2637/100', 1, 25, '1')],
        },
        TAIL={
            2: ('2637/100', '26', '7/4', '26'),
            3: ('367/25', '1467/100', '281/100', '339/10'),
            4: ('3047/100', '152/5', '14/5', '70'),
            5: ('367/25', '73/5', '13/4', '40'),
            6: ('2637/100', '26', '2', '30'),
        },
        LOWER='lower_tests_n4.json',
        LOWER_BREAKS=[0, 1, 2, 4, 8, 16, 32, 64, 112, 160],
        LOWER_DEG=12,
    ),
}


def asymptotic_coefficients(n):
    """Coefficients (a_1, a_2, a_3) in f_n = 1 - a_1 r^-2 - a_2 r^-4 - a_3 r^-6 + O(r^-8),
    as stated in del Pino-Liu-Musso-Wei-Yang, arXiv:2608.04976v1, eq. (asym)."""
    from fractions import Fraction as F
    return (F(n * n, 2), F(n * n * (n * n + 8), 8), F(n * n * (n ** 4 + 32 * n * n + 128), 16))
