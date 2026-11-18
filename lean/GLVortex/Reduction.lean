import GLVortex.Profile
import GLVortex.Gram

set_option autoImplicit false

namespace GLVortex

/-- Channels indexed 0,...,4 represent m=2,...,6. -/
def K (i : Fin 5) : ℕ := if i.val = 2 then 2 else 1

/-- The two-dimensional negative subspace in Lemma 5.
Negativity of the two diagonal entries alone would not suffice. -/
theorem negative_definite_two_by_two (a b c : ℝ)
    (ha : a < 0) (hdet : 0 < a*c-b^2)
    (u v : ℝ) (hne : u ≠ 0 ∨ v ≠ 0) :
    a*u^2 + 2*b*u*v + c*v^2 < 0 := by
  have hid : a*(a*u^2+2*b*u*v+c*v^2) =
      (a*u+b*v)^2 + (a*c-b^2)*v^2 := by ring
  by_cases hv : v = 0
  · have hu : u ≠ 0 := hne.resolve_right (by simpa [hv])
    have hu2 : 0 < u^2 := sq_pos_of_ne_zero hu
    simp only [hv, zero_pow (by decide : 2 ≠ 0), mul_zero, add_zero]
    exact mul_neg_of_neg_of_pos ha hu2
  · have hv2 : 0 < v^2 := sq_pos_of_ne_zero hv
    have hp : 0 < (a*c-b^2)*v^2 := mul_pos hdet hv2
    have hs := sq_nonneg (a*u+b*v)
    have hprod : 0 < a*(a*u^2+2*b*u*v+c*v^2) := by linarith
    nlinarith

/-- The exact Gram entries (not the rounded decimal table) define a negative form. -/
theorem lemma_5_central_gram (u v : ℝ) (hne : u ≠ 0 ∨ v ≠ 0) :
    (gram4_00 : ℝ)*u^2 + 2*(gram4_01 : ℝ)*u*v + (gram4_11 : ℝ)*v^2 < 0 := by
  apply negative_definite_two_by_two _ _ _ _ _ u v hne
  · exact_mod_cast gram4_negative
  · exact_mod_cast gram4_det_positive

/-- Lemmas 4 and 5 sandwich the index. The `nunull` identity is an
explicit upstream hypothesis; no spectral theory is hidden in the proof. -/
theorem channel_sandwich (index nullity nu bound : ℕ)
    (nunull : nu = index + nullity)
    (trace_upper : nu ≤ bound) (negative_lower : bound ≤ index) :
    index = bound ∧ nullity = 0 := by omega

/-- All finite upper sign conditions used in Lemma 4. Identities relating these
polynomials to derivatives/potentials and the trace theorem are external inputs. -/
def UpperConditions : Prop :=
  (∀ x, 0 ≤ x → 0 < denominatorL x) ∧
  (∀ x, 0 ≤ x → 0 < channelA2 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA3 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA4 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA5 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA6 x) ∧
  (∀ x, 0 ≤ x → x ≤ 2637/100 → 0 < inner2_0 x ∧ 0 < innerResidual2_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ 367/25 → 0 < inner3_0 x ∧ 0 < innerResidual3_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ 144579/25000 → 0 < inner4_0 x ∧ 0 < innerResidual4_0 x) ∧
  (∀ x, 144579/25000 ≤ x → x ≤ 3047/100 → 0 < inner4_1 x ∧ 0 < innerResidual4_1 x) ∧
  (∀ x, 0 ≤ x → x ≤ 367/25 → 0 < inner5_0 x ∧ 0 < innerResidual5_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ 2637/100 → 0 < inner6_0 x ∧ 0 < innerResidual6_0 x) ∧
  (∀ x, 2637/100 ≤ x → 0 < tail2 x) ∧
  (∀ x, 367/25 ≤ x → 0 < tail3 x) ∧
  (∀ x, 3047/100 ≤ x → 0 < tail4 x) ∧
  (∀ x, 367/25 ≤ x → 0 < tail5 x) ∧
  (∀ x, 2637/100 ≤ x → 0 < tail6 x)

theorem upper_conditions : UpperConditions := by
  refine ⟨denominatorL_pos, channelA2_pos, channelA3_pos, channelA4_pos,
    channelA5_pos, channelA6_pos, ?_, ?_, ?_, ?_, ?_, ?_,
    tail2_pos, tail3_pos, tail4_pos, tail5_pos, tail6_pos⟩
  all_goals intro x ha hb
  · exact ⟨inner2_0_pos x ha hb, innerResidual2_0_pos x ha hb⟩
  · exact ⟨inner3_0_pos x ha hb, innerResidual3_0_pos x ha hb⟩
  · exact ⟨inner4_0_pos x ha hb, innerResidual4_0_pos x ha hb⟩
  · exact ⟨inner4_1_pos x ha hb, innerResidual4_1_pos x ha hb⟩
  · exact ⟨inner5_0_pos x ha hb, innerResidual5_0_pos x ha hb⟩
  · exact ⟨inner6_0_pos x ha hb, innerResidual6_0_pos x ha hb⟩

/-- Analytic hypotheses for Theorem 1.
The kernel propositions must be instantiated with the actual PDE kernel
equalities; no PDE operator or Fourier expansion is constructed in Lean. -/
structure AnalyticInputs (f : ℝ → ℝ) (boundedKernelSymmetries fullL2Trivial : Prop)
    (index nullity nu : Fin 5 → ℕ) (fullIndex : ℕ) where
  profile : ProfileFacts f
  nunull : ∀ i, nu i = index i + nullity i
  /-- Ratio comparison, monotonicity, scalar reduction, interval tiling and
  trace-point proposition; derivative identities are not formalized here. -/
  trace_upper : ProfileFacts f → PolynomialConditions → UpperConditions → ∀ i, nu i ≤ K i
  /-- Majorant/form-integral identification, admissibility, monotonicity and
  variational principle. Exact Gram signs are arguments, not assumed signs. -/
  negative_lower : ProfileFacts f → PolynomialConditions →
    gram2_00 < 0 → gram3_00 < 0 → gram4_00 < 0 →
    0 < gram4_00*gram4_11-gram4_01^2 → gram5_00 < 0 → gram6_00 < 0 →
    ∀ i, K i ≤ index i
  /-- General-degree bounded-to-L² transfer plus phase/translation/high-channel
  classification and Fourier reconstruction. -/
  bounded_transfer : (∀ i, nullity i = 0) → boundedKernelSymmetries
  /-- Symmetry fields are not L²; other channels have trivial L² kernels. -/
  l2_reconstruction : (∀ i, nullity i = 0) → fullL2Trivial
  /-- Two real copies; all remaining channels have zero index. -/
  angular_multiplicity : fullIndex = 2 * ∑ i, index i

/-- Lemma 4 assumes the analytic trace reduction. -/
theorem lemma_4_upper {f : ℝ → ℝ} {B L : Prop} {ind nul nu : Fin 5 → ℕ} {full : ℕ}
    (h : AnalyticInputs f B L ind nul nu full) : ∀ i, nu i ≤ K i :=
  h.trace_upper h.profile all_polynomial_conditions upper_conditions

/-- Lemma 5, conditional on admissibility and exact identification of the saved
rational form integrals with the variational test subspaces. -/
theorem lemma_5_lower {f : ℝ → ℝ} {B L : Prop} {ind nul nu : Fin 5 → ℕ} {full : ℕ}
    (h : AnalyticInputs f B L ind nul nu full) : ∀ i, K i ≤ ind i :=
  h.negative_lower h.profile all_polynomial_conditions gram2_negative gram3_negative gram4_negative
    gram4_det_positive gram5_negative gram6_negative

/-- Theorem 1: exact five-channel spectral sandwich, trivial intermediate
kernels, conditional full-kernel conclusions and real Morse index 12. -/
theorem theorem_1 {f : ℝ → ℝ} {B L : Prop} {ind nul nu : Fin 5 → ℕ} {full : ℕ}
    (h : AnalyticInputs f B L ind nul nu full) :
    (∀ i, ind i = K i ∧ nul i = 0) ∧ B ∧ L ∧ full = 12 := by
  have hc : ∀ i, ind i = K i ∧ nul i = 0 := by
    intro i
    exact channel_sandwich _ _ _ _ (h.nunull i) (lemma_4_upper h i) (lemma_5_lower h i)
  refine ⟨hc, h.bounded_transfer (fun i => (hc i).2),
    h.l2_reconstruction (fun i => (hc i).2), ?_⟩
  rw [h.angular_multiplicity]
  have heq : ind = K := funext (fun i => (hc i).1)
  rw [heq]
  norm_num [K, Fin.sum_univ_succ]

end GLVortex
