import Mathlib.Tactic

set_option autoImplicit false

namespace GLVortex

/-- Homogeneous power/Bernstein evaluation; binomial weights are included in `cs`. -/
def heval : List ℝ → ℝ → ℝ → ℝ
  | [], _, _ => 0
  | c :: cs, u, v => c * v ^ cs.length + u * heval cs u v

/-- A nonempty strictly positive coefficient list gives a strictly positive
homogeneous polynomial whenever both variables are nonnegative and not both zero. -/
theorem heval_pos (cs : List ℝ) (hne : cs ≠ [])
    (hc : ∀ c ∈ cs, 0 < c) (u v : ℝ)
    (hu : 0 ≤ u) (hv : 0 ≤ v) (huv : 0 < u + v) : 0 < heval cs u v := by
  induction cs with
  | nil => exact False.elim (hne rfl)
  | cons c cs ih =>
    have hc0 : 0 < c := hc c (by simp)
    have htail : ∀ d ∈ cs, 0 < d := by
      intro d hd
      exact hc d (by simp [hd])
    cases cs with
    | nil => simpa [heval] using hc0
    | cons d ds =>
      have hp := ih (by simp) htail
      have hnon : 0 ≤ heval (d :: ds) u v := le_of_lt hp
      by_cases hvp : 0 < v
      · have hfirst : 0 < c * v ^ (d :: ds).length :=
          mul_pos hc0 (pow_pos hvp _)
        exact add_pos_of_pos_of_nonneg hfirst (mul_nonneg hu hnon)
      · have hv0 : v = 0 := le_antisymm (le_of_not_gt hvp) hv
        have hup : 0 < u := by linarith
        have hsecond : 0 < u * heval (d :: ds) u v := mul_pos hup hp
        have hfirst : 0 ≤ c * v ^ (d :: ds).length :=
          mul_nonneg (le_of_lt hc0) (pow_nonneg hv _)
        exact add_pos_of_nonneg_of_pos hfirst hsecond

end GLVortex
