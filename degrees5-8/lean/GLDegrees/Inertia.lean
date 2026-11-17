import GLDegrees.Finite
namespace GLDegrees

theorem inertiaCount_append (xs ys : List (List ℤ)) :
    inertiaCount (xs ++ ys) = (do
      let a ← inertiaCount xs
      let b ← inertiaCount ys
      pure (a+b)) := by
  induction xs with
  | nil => cases h : inertiaCount ys <;> simp [inertiaCount, h]
  | cons x xs ih =>
    simp only [List.cons_append, inertiaCount, ih]
    cases hx : boxNegative x <;> cases hs : inertiaCount xs <;>
      cases ht : inertiaCount ys <;> simp [hx, hs, ht, Nat.add_assoc]

#print axioms inertiaCount_append
end GLDegrees
