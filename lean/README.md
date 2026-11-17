# Degree-four vortex: finite algebra and conditional spectral reduction

Lean 4.35.0; Mathlib v4.35.0. Build with `lake build`.

## PROVED in Lean

- `heval_pos`: positivity of homogeneous power/Bernstein evaluation with strictly positive coefficients on the nonnegative quadrant excluding the origin.
- `Certificates.lean`: 40 real-polynomial positivity theorems. These cover both barrier denominators, both residual signs, all five A and E denominators, all six inner polynomials and their normalized residuals, all five tail polynomials, and all nine lower-majorant pieces. Each proof checks a polynomial identity with `ring` and exact rational coefficient inequalities with `norm_num`; closed finite intervals and the entire half-line are covered, not samples.
- `all_polynomial_conditions`: the conjunction of those 40 conclusions.
- `Gram.lean`: signs of the unrounded exact-rational negative-form entries and the central 2-by-2 determinant. These are the exact values extracted from the saved test functions, not the decimal table rounded back to rational numbers.
- `negative_definite_two_by_two` and `lemma_5_central_gram`: determinant criterion implies strict negativity for every nonzero real pair, not merely for the two coordinate vectors.
- `Algebra.lean`: the nonnegative monotonicity factor, exact scalar square completion, and its pointwise lower-bound consequence.
- `channel_sandwich`: upper nonpositive dimension and lower negative dimension plus the explicit index/nullity identity give exact index and zero nullity.
- `theorem_1`: conditional five-channel index tuple (1,1,2,1,1), vanishing intermediate nullities, conditional full kernel conclusions, and the real index 12.
- Numbered interfaces `lemma_2_profile_enclosure`, `lemma_4_upper`, `lemma_5_lower` mirror manuscript Lemmas 2, 4 and 5. Proposition 3's analytic trace theorem is an input, not a Lean theorem.

## Explicit hypotheses, NOT proved by Lean

`ProfileFacts` states the classical profile ODE, differentiability, positivity, origin expansion with positive leading coefficient, and the infinity expansion with a rigorous O(r^-8) remainder. No profile existence, uniqueness, expansion remainder or shooting accuracy is asserted by Lean.

`AnalyticInputs` exposes the upstream index/nullity identity; ratio/profile comparison, form monotonicity, scalar reduction and the trace-point theorem; identification of polynomial residuals with the differential expressions; interval tiling in logarithmic radius; form-domain admissibility; identification of the extracted exact Gram entries with the integrals of the specified test functions; the variational principle; bounded-to-L2 transfer; phase/translation/high-channel classification; Fourier reconstruction; and the two real angular copies. The full bounded-kernel and L2-kernel propositions are abstract parameters to be instantiated by the PDE statements. **This is not a formal construction of the PDE operator or an unconditional Lean proof of vortex nondegeneracy.**

The saved polynomial and Gram constants are proposals checked by Lean. Their extraction from analytic differential formulas and trial integrals remains part of the external certificate verification. Agreement of generated data and a kernel-checked sign does not formalise that extraction.

## Reproduction and axioms

The checked-in `.lean` files are sufficient; no Python generator is needed. Both Lean and Mathlib are pinned to v4.35.0. From this directory:

```
lake update
lake exe cache get
lake build
lake env lean AxiomCheck.lean
```

`GLVortex.lean` prints the axioms of every main theorem and the aggregate polynomial certificate. The proofs use `norm_num`, `ring`, linear arithmetic and kernel proof terms. They do not use `native_decide`, custom axioms or proof holes. Standard Lean/Mathlib axioms may appear in the printed dependency lists.

The project uses a portable git dependency pinned to Mathlib v4.35.0. Fetch prebuilt Mathlib oleans with `lake exe cache get` rather than rebuilding Mathlib from source.
The portable Git-based dependency configuration and CI have not been tested in a clean clone. The checked results used prebuilt dependencies; these commands are reproduction instructions, not evidence of an additional portable build.

## Higher-degree boundary

This is the checked degree-four project only. `theorem_1` is conditional bookkeeping, not an unconditional PDE theorem. In the separate `../degrees5-8/lean/` project, 15 of 97 files (96 generated modules plus the root import file) are checked; the remainder are not machine-checked. Its root imports unchecked modules and does not build as shipped. `Tail5M2` was accepted only with `-M 10000`, not the earlier `-M 5000` limit; the shipped higher-degree lakefile now sets `-M 10000`. The root README gives the exact invocation. No complete higher-degree certificate-chain theorem is claimed. The root README lists the checked modules.

