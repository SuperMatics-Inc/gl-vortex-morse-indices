# Degree-four vortex: finite algebra and conditional spectral reduction

Lean 4.35.0; Mathlib v4.35.0. Build with `lake build`.

## Results checked in Lean

- `heval_pos`: positivity of homogeneous power/Bernstein evaluation with strictly positive coefficients on the nonnegative quadrant excluding the origin.
- `Certificates.lean`: 40 real-polynomial positivity theorems. These cover both barrier denominators, both residual signs, all five A and E denominators, all six inner polynomials and their normalized residuals, all five tail polynomials, and all nine lower-majorant pieces. Each proof checks a polynomial identity with `ring` and exact rational coefficient inequalities with `norm_num`; closed finite intervals and the entire half-line are covered, not samples.
- `all_polynomial_conditions`: the conjunction of those 40 conclusions.
- `Gram.lean`: signs of the unrounded exact-rational negative-form entries and the central 2-by-2 determinant. These are the exact values extracted from the saved test functions, not the decimal table rounded back to rational numbers.
- `negative_definite_two_by_two` and `lemma_5_central_gram`: determinant criterion implies strict negativity for every nonzero real pair, not merely for the two coordinate vectors.
- `Algebra.lean`: the nonnegative monotonicity factor, exact scalar square completion, and its pointwise lower-bound consequence.
- `channel_sandwich`: upper nonpositive dimension and lower negative dimension plus the explicit index/nullity identity give exact index and zero nullity.
- `theorem_1`: conditional five-channel index tuple (1,1,2,1,1), vanishing intermediate nullities, conditional full kernel conclusions, and the real index 12.
- The statements `lemma_2_profile_enclosure`, `lemma_4_upper`, `lemma_5_lower` mirror manuscript Lemmas 2, 4 and 5. Proposition 3's analytic trace theorem is an input, not a Lean theorem.

## Analytic hypotheses

`ProfileFacts` assumes the classical profile ODE, differentiability, positivity and the origin expansion with positive leading coefficient. It also assumes the expansion at infinity with a rigorous O(r^-8) remainder. Existence, uniqueness and shooting accuracy lie outside the formalization.

`AnalyticInputs` assumes the index/nullity identity, profile ratio comparison, form monotonicity, scalar reduction and the trace-point theorem. It also assumes the identification of polynomial residuals with differential expressions, the logarithmic interval tiling and form-domain admissibility. The Gram integral identities and variational principle are inputs. So are the bounded-to-L2 transfer, the classification of phase, translation and high-channel modes, Fourier reconstruction and the two real angular copies. The bounded-kernel and L2-kernel statements are parameters for the corresponding PDE results.

Lean checks signs of the supplied polynomial and Gram constants. Their extraction from differential formulas and trial integrals is part of the external certificate verification.

## Reproduction and axioms

The checked-in `.lean` files suffice; no Python generator is needed. Both Lean and Mathlib use version 4.35.0. From this directory:

```
lake update
lake exe cache get
lake build
lake env lean AxiomCheck.lean
```

`GLVortex.lean` prints the axioms of every main theorem and the aggregate polynomial certificate. The proofs use `norm_num`, `ring`, linear arithmetic and kernel proof terms. They do not use `native_decide`, custom axioms or proof holes. Standard Lean/Mathlib axioms may appear in the printed dependency lists.

Fetch Mathlib’s prebuilt `.olean` files with `lake exe cache get`. The checks used prebuilt dependencies.

## Degrees five through eight

The separate `../degrees5-8/lean/` project has checks for 15 of its 97 files (96 generated modules and the root import file). The root imports unchecked modules and does not build as supplied. `Tail5M2` passed with `-M 10000` after exceeding the earlier `-M 5000` limit. The higher-degree lakefile uses `-M 10000`. The main README lists the checked modules and gives the exact invocation.

