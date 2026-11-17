import GLDegrees.Finite
import Mathlib.Algebra.Polynomial.CoeffList
set_option autoImplicit false
namespace GLDegrees

theorem eqP_sound {p q : IP} {n : ℕ} (h : eqP p q n = true) :
    Polynomial.ofCoeffList p = Polynomial.ofCoeffList q := by
  simp only [eqP, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true, beq_iff_eq] at h
  ext i
  simp only [Polynomial.coeff_ofCoeffList, List.getD_eq_getElem?_getD]
  by_cases hi : i < n
  · exact h.2 i (List.mem_range.mpr hi)
  · have hp : p[i]? = none := List.getElem?_eq_none_iff.mpr (Nat.le_trans h.1.1 (Nat.le_of_not_gt hi))
    have hq : q[i]? = none := List.getElem?_eq_none_iff.mpr (Nat.le_trans h.1.2 (Nat.le_of_not_gt hi))
    simp [hp, hq]

#print axioms eqP_sound
end GLDegrees
