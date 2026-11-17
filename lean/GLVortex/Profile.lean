import GLVortex.FiniteConditions
import Mathlib.Analysis.Calculus.Deriv.Basic

namespace GLVortex

/-- Classical profile facts are hypotheses, not consequences of finite algebra.
The two quantified bounds express the rigorous O(r^6) origin and O(r^-8)
infinity remainders; no numerical shooting data enters this definition. -/
structure ProfileFacts (f : ℝ → ℝ) : Prop where
  positive : ∀ r, 0 < r → 0 < f r
  smooth : ∀ r, 0 < r → DifferentiableAt ℝ f r ∧ DifferentiableAt ℝ (deriv f) r
  ode : ∀ r, 0 < r →
    deriv (deriv f) r + deriv f r / r - 16 * f r / r^2 + (1-(f r)^2)*f r = 0
  origin : ∃ κ C δ : ℝ, 0 < κ ∧ 0 ≤ C ∧ 0 < δ ∧
    ∀ r, 0 < r → r ≤ δ → |f r - κ*r^4| ≤ C*r^6
  infinity_remainder : ∃ C R : ℝ, 0 ≤ C ∧ 0 < R ∧
    ∀ r, R ≤ r → |f r - (1 - 8/r^2 - 48/r^4 - 896/r^6)| ≤ C/r^8

def ProfileEnclosure (f : ℝ → ℝ) : Prop :=
  ∀ r, 0 < r → r^8 / denominatorL (r^2) < (f r)^2 ∧
    (f r)^2 < r^8 / denominatorU (r^2)

/-- Lemma 2. Ratio comparison (del Pino et al., Lemma 2.2) and the
identification of the residual polynomials with E_4(g) remain hypotheses. -/
theorem lemma_2_profile_enclosure (f : ℝ → ℝ) (hf : ProfileFacts f)
    (ratio_comparison : ProfileFacts f → PolynomialConditions → ProfileEnclosure f) :
    ProfileEnclosure f := ratio_comparison hf all_polynomial_conditions

end GLVortex
