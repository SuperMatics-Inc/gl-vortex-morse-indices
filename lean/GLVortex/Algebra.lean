import GLVortex.Positivity

namespace GLVortex

/-- Pointwise nonnegative factor in the form monotonicity identity. -/
theorem monotonicity_factor (a b : ℝ) : 0 ≤ a^2+a*b+b^2 := by
  nlinarith [sq_nonneg (a+b), sq_nonneg (a-b)]

/-- The exact square completion in equation (15), independent of all ODE facts. -/
theorem scalar_square_completion (A W c d : ℝ) (hne : A+W ≠ 0) :
    A*c^2+W*(c-d)^2 =
      (A+W)*(c-W*d/(A+W))^2 + A*W*d^2/(A+W) := by
  field_simp [hne]
  ring

/-- The algebraic scalar lower bound; integration and domain justification
are separate analytic inputs. -/
theorem scalar_lower_bound (A W c d : ℝ) (hA : 0 < A) (hW : 0 ≤ W) :
    A*W*d^2/(A+W) ≤ A*c^2+W*(c-d)^2 := by
  have hsum : 0 < A+W := by linarith
  rw [scalar_square_completion A W c d (ne_of_gt hsum)]
  have hs : 0 ≤ (A+W)*(c-W*d/(A+W))^2 := mul_nonneg (le_of_lt hsum) (sq_nonneg _)
  linarith

end GLVortex
