import Mathlib.Tactic
set_option autoImplicit false
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace GLDegrees
abbrev IP := List ℤ
def coeff (p : IP) (i : ℕ) : ℤ := (p[i]?).getD 0
def add (p q : IP) : IP := (List.range (max p.length q.length)).map (fun i => coeff p i + coeff q i)
def scale (c : ℤ) (p : IP) : IP := p.map (fun a => c*a)
def mul (p q : IP) : IP :=
  (List.range (p.length+q.length)).map (fun k =>
    ((List.range (k+1)).map (fun i => coeff p i * coeff q (k-i))).sum)
def deriv (p : IP) : IP := (List.range p.length).map (fun (i : ℕ) => ((i+1 : ℕ) : ℤ)*coeff p (i+1))
def shift (k : ℕ) (p : IP) : IP := List.replicate k 0 ++ p
def eqP (p q : IP) (n : ℕ) : Bool :=
  decide (p.length ≤ n ∧ q.length ≤ n) && (List.range n).all (fun i => coeff p i == coeff q i)
def mulLinear (a b : ℤ) (p : IP) : IP :=
  (List.range (p.length+1)).map (fun i => a*coeff p i + if i=0 then 0 else b*coeff p (i-1))
def powLinear (a b : ℤ) (n : ℕ) : IP :=
  (List.range (n+1)).map (fun i => (n.choose i : ℤ)*a^(n-i)*b^i)
def homogeneous (u0 u1 v0 v1 : ℤ) : IP → IP
  | [] => []
  | c::cs => add (scale c (powLinear v0 v1 cs.length)) (mulLinear u0 u1 (homogeneous u0 u1 v0 v1 cs))
def positiveCoefficients (p : IP) : Bool := !p.isEmpty && p.all (fun c => decide (0<c))
def residual (n : ℕ) (den : ℤ) (h : IP) : IP :=
  let z : ℤ := n
  let dh := deriv h
  let b := add (add (scale (-z) (mul [2*(z+1),-(z+2)] h))
      (scale 4 (mul (mul [1,-1] [z+1,-(z+2)]) dh)))
      (scale 4 (mul [0,1,-2,1] (deriv dh)))
  add (add (scale (z*z*den*den) h) (scale (-(z*z)) (shift n (mul (mul h h) h))))
    (scale (den*den) (mul [1,-1] b))
def gramDet3 (a b c d e f : ℤ) : ℤ := a*d*f+2*b*c*e-a*e*e-d*c*c-f*b*b
def composeScaled (a b den : ℤ) : IP → IP
  | [] => []
  | c::cs => add [c*den^cs.length] (mulLinear a b (composeScaled a b den cs))
def tailW (n : ℕ) (h : IP) : IP := scale ((n : ℤ)^2) (shift (n+1) (mul h h))
def tailA (n : ℕ) (k den : ℤ) (h : IP) : IP :=
  add [k*k*den^2, -(k*k+(n : ℤ)^2)*den^2] (scale 2 (tailW n h))
def trim (p : IP) : IP := (p.reverse.dropWhile (fun c => c == 0)).reverse
def fourMin (a b c d : ℤ) : ℤ := min (min a b) (min c d)
def fourMax (a b c d : ℤ) : ℤ := max (max a b) (max c d)
def boxNegative : List ℤ → Option ℕ
  | [al,ah,bl,bh,cl,ch] =>
    if al≤ah ∧ bl≤bh ∧ cl≤ch then
      let dl := fourMin (al*cl) (al*ch) (ah*cl) (ah*ch) - fourMax (bl*bl) (bl*bh) (bh*bl) (bh*bh)
      let du := fourMax (al*cl) (al*ch) (ah*cl) (ah*ch) - fourMin (bl*bl) (bl*bh) (bh*bl) (bh*bh)
      if du<0 then some 1 else if 0<dl ∧ ah<0 then some 2 else if 0<dl ∧ 0<al then some 0 else none
    else none
  | _ => none
def inertiaCount : List (List ℤ) → Option ℕ
  | [] => some 0
  | b::bs => do
    let c ← boxNegative b
    let d ← inertiaCount bs
    pure (c+d)
end GLDegrees
