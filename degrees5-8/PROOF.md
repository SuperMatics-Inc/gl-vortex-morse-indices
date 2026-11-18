# Variational bounds for the vortex index

Sections 2–5 give the argument for the passage from the finite certificates to the differential operator, under the analytic hypotheses in Section 1. This passage is a written proof rather than part of the Lean formalization.

## 1. Imported analytic hypotheses

The analytic reference throughout is del Pino–Liu–Musso–Wei–Yang,
arXiv:2608.04976v1. Its Theorems 1.1–1.2
compute degrees 2 and 3 only. Their conclusions are not imported for degrees 5--8.
Our explicit analytic hypotheses are the standard positive, strictly monotone,
regular profile, its origin expansion and limit one as in (1.1)--(1.2);
the real Fourier decomposition (2.1), with two copies of every m>=1;
the closed-form/min--max setting of Lemma 2.1 and (2.2);
and the ratio-comparison principle Lemma 2.2, **including its ordering near infinity**.
We use the general-n identities (3.1), (3.2)--(3.6), and (3.7) for low and high
channels, not the degree-specific numerical conclusions of that paper.
Their proofs are degree-independent: (3.1) expands a square using the profile
equation; q=f'+nf/r is positive, and p=nf/r-f'>0 follows because
z=rf'/f starts below n and satisfies dz/d(log r)=n²-r²(1-f²)-z²,
so its first contact with n would have strictly negative derivative.
Differentiating the vortex gives the positive pair (q,p) in (3.4)--(3.5).
Expanding two squares gives (3.6), hence Q1>=0. Subtracting angular squares gives
Qm-Q1=integral[(m-1)(m-2n+1)a²+(m-1)(m+2n+1)b²]dr/r.
This is nonnegative for every m>=2n-1. Q0>=0 follows from (3.1).
These arguments apply to every positive integer n; no degree-2/3 barrier is used.
Extension by form closure, regular-origin traces and the variational interpretation
remain analytic hypotheses, not Lean theorems.

For intermediate m=2,…,2n−2 put k=|n−m|, j=n+m and t=log r, x=e^(2t). The channel form is

`Q_f(v)=∫ (|v_t|² + vᵀ B_f(t) v) dt`,
`B_f=diag(k²,j²)−x I+x f² [[2,1],[1,2]]`.

The last constant matrix has eigenvalues1,3. Consequently `g_L≤f≤g_U` implies `Q_gL≤Q_f≤Q_gU`. No kinetic energy of the massive component is discarded, unlike the inherited scalar Schur bound.

## 2. Exact compactified profile barriers

Set s=n², y=x/(s+x) and `g=y^(n/2) H(y)`, with rational polynomial H. All stored H satisfy `H(1)=1`, `H'(1)=−(n−1)/2`. If E denotes the profile differential operator, direct differentiation gives

`s E(g)/y^(n/2) = R(H)`,
`R(H)=s H(1−y^n H²)+(1−y)[−n(2(n+1)−(n+2)y)H +4(1−y)((n+1)−(n+2)y)H' +4y(1−y)²H'']`.

The finite checks establish H>0 and `R(H_L)/(1−y)²>0`, `−R(H_U)/(1−y)²>0` on [0,1]. They also check `nH+2yH'>0` (profile monotonicity), and `(H_U−H_L)/(1−y)²>0`. Numerical BVP digits only design H; they are not trusted inputs to acceptance.

### Infinity ordering and a rigorous remainder

We first establish the ordering at infinity required by Lemma 2.2. Write s=n², u=r^-2. Direct Taylor expansion gives
`g²=1−s u+c_g u²+O(u³)`,
`c_g=s²[H''(1)+(5−n²)/4]`.
For the actual profile,
`f=1−s/(2r²)−s(s+8)/(8r⁴)−s(s²+32s+128)/(16r⁶)+O(r^-8)`, hence
`f²=1−s u−2s u²+O(u³)`.
To bound the remainder, let F=1+a r^-2+b r^-4+d r^-6 with a=−s/2, b=−s(s+8)/8 and d=−s(s²+32s+128)/16.
Substitution into E shows the r^-2, r^-4 and r^-6 coefficients vanish, so E(F)=O(r^-8).
For w=f−F, exact subtraction gives
`(-d²/dr²−r^-1 d/dr+c(r))w=E(F)`,
`c=s/r²+f²+fF+F²−1 -> 2`.
Choose R so c>=1 and R²>=128. For h=C r^-8,
`(-d²/dr²−r^-1 d/dr+c)h=C(c−64/r²)r^-8>=C r^-8/2`.
Choose C large enough to dominate |E(F)| and |w(R)|R^8. Since w,h tend to zero,
the exterior maximum principle applied to w−h and −w−h gives |w|<=h:
any positive maximum is attained at a finite interior point, where the
differential inequality contradicts c>0. This proves the O(r^-8)
remainder using only the profile equation and f->1. At n=4 this gives the
coefficient −896 needed for the degree-four barriers.

Let red be the **signed positive** polynomial stored by the certificate.
The exact identity
`s red(1)=c_g+2s` for the upper barrier and
`s red(1)=−(c_g+2s)` for the lower barrier
is reconstructed and checked by `checked_barrier`. Thus c_L<−2s<c_U strictly.
The expansions give g_L<f<g_U near infinity. At the origin,
g is a positive constant times r^n(1+O(r²)), since H(0)>0.
All hypotheses of the imported ratio-comparison Lemma 2.2 now hold, yielding
the global enclosure.

## 3. Exact lower-index witnesses

Partition a finite t interval using exact rational x endpoints. Its lengths obey

`ell = atanh(z)`, `z=(b−a)/(b+a)`,
`z+z³/3+z⁵/5 ≤ ell ≤ z+z³/3+z⁵/5+z⁷/[7(1−z²)]`.

The code rounds these bounds outwards to denominator 10^15. Since x and g_U² are increasing,
`B_gU(t) ≤ Bplus=diag(k²,j²)−aI+b g_U(b)² M`.
For vector linear interpolation L between nodal values, the kinetic contribution is bounded above using 1/ell_L. The potential integral is bounded above by using ell_U and adding `(ell_U−ell_L)||Bplus||/ell_U` to its diagonal. The exact polynomial interpolation integrals have diagonal weight 1/3 and off-diagonal weight 1/6.

The resulting symmetric 4×4 cell is rounded to denominator 10^12. Each entry changes by less than1/10^12; adding4/10^12 to its diagonal gives a Loewner upper bound by the row-sum bound. The right endpoint is zero, giving compact support. At the left endpoint x0, extend each component by `(x/x0)^(k/2)` and `(x/x0)^(j/2)`. Its kinetic/angular contribution is exactly diag(k,j). Dropping the negative potential and using `M≤3I` bounds the remaining origin contribution by
`diag(3x0 g_U(x0)²/[2(k+1)],3x0 g_U(x0)²/[2(j+1)])`.

Stored nodal vectors have common denominator 10^7. Their upper Gram matrix is therefore integer-valued after multiplication by 10^26. Strict alternating signs of all leading principal minors certify a negative-definite trial subspace by Sylvester's criterion. This certifies an index lower bound only, without needing a numerical eigenvalue enclosure.

## 4. Coupled upper-index witnesses

### Origin and tail

On x≤x0≤1, ignore the positive profile matrix. For scalar angular degree q≥0 use `h=x^(q/2) exp(−x/[2(q+1)])`. Its t-log derivative is `s=q−x/(q+1)` and
`(q²−x)−h_tt/h = x−x²/(q+1)² ≥0`.
The ground-state identity gives the origin lower form `s(x0)|v(t0)|²`. At q=0 the boundary term at−∞ vanishes despite h approaching1; this requires the standard finite-energy origin regularity hypothesis.

For the tail x≥X, choose an integer L≥0. The verifier proves `B_gL−L² I` pointwise positive definite: after multiplication by1−y, its first diagonal is positive and determinant divided by1−y is positive on[X/(n²+X),1]. Completing a square on the half-line gives
`∫[|v'|²+vᵀB_gL v]dt ≥ ∫[|v'|²+L²|v|²]dt = L|v(tX)|²+∫|v'+Lv|²dt ≥ L|v(tX)|²`.
For L>0, the boundary term at infinity vanishes by H¹∩L² tail regularity. For L=0, the contribution is nonnegative. Add `L I` to the last nodal block, retaining the free endpoint at X. The certificates choose L=floor(m/2) and check both polynomial inequalities for every channel. Any nonnegative L satisfying the same inequalities works. The earlier choice L=0 at R20 did not resolve n=7, m=6.

### Interior bubbles

Write v=L+z on each cell, where z vanishes at both endpoints. The derivative cross term integrates to zero. For a lower potential D(t) with `||D||≤C`, Young's inequality and `∫|z|²≤ell²/π²∫|z'|²≤ell_U²/9∫|z'|²` give

`Q_D(L+z) ≥ ∫|L'|² + ∫Lᵀ D L − eta∫|L|²`
provided `eta≥C² ell_U²/(9−C ell_U²)>0`.

For the endpoint-constant scheme, take `D=diag(k²,j²)−bI+a g_L(a)² M` and integrate exactly.

Chord scheme: interpolate B_gL linearly between endpoint matrices, with endpoint xg_L² rounded down. Put G=y^n H_L². Direct differentiation gives
`(xG)_tt=4x Z`,
`Z=G+y(1−y)(3−2y)G'+y²(1−y)²G''`.
The exact polynomial tests `10±Z>0` imply `||B_gL''||≤4b(1+30)` on the cell. Subtract `delta=ell_U² b·31/2` times I from the chord. The standard second-derivative interpolation bound makes this D a pointwise lower bound. Its norm is bounded by the maximum endpoint row-sum norm plus delta.

Integrating Lᵀ times the linear chord gives endpoint diagonal weights(1/4,1/12) and off-diagonal weights(1/12,1/12). Subtracting delta+eta uses weights 1/3 and 1/6. A further outward diagonal correction handles ell_L≤ell≤ell_U; kinetic energy is bounded below using1/ell_U. These are all exact rational finite matrices.

### Inertia

Assemble the 2×2 block tridiagonal matrix including the origin Robin term. Every block Schur pivot is enclosed by rational intervals rounded outwards to the certificate's fixed denominator. The independent verifier checks each enclosure against the previous enclosure and the exact cell matrices. Nonzero determinant intervals and the first diagonal sign give the exact pivot inertia; congruence then gives the whole finite-matrix inertia.

Any negative subspace of the continuous form injects into the nodal values: the discarded bubble form is nonnegative, and a vector with zero nodal values cannot be strictly negative. Restricting to the inverse image of the finite matrix's nonnegative spectral subspace proves `index(Q_f)≤index(matrix)`. Thus matching exact upper and lower channel counts certify the full negative index, conditional on §1 and this written variational bridge. No kernel classification or positive full-plane gap is claimed.

## 5. Scope

Floating-point eigenvalues on a finite box are numerical estimates. The full-plane essential spectrum begins at zero. In particular, the third level at n=7, m=6 lies near threshold in the finite-box calculation. The certificates count negative directions and give no positive full-plane gap.
