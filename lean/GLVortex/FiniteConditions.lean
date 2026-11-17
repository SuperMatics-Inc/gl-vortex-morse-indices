import GLVortex.Certificates

namespace GLVortex

def PolynomialConditions : Prop :=
  (∀ x, 0 ≤ x → 0 < denominatorU x) ∧
  (∀ x, 0 ≤ x → 0 < denominatorL x) ∧
  (∀ x, 0 ≤ x → 0 < negativeResidualU x) ∧
  (∀ x, 0 ≤ x → 0 < residualL x) ∧
  (∀ x, 0 ≤ x → 0 < channelA2 x) ∧
  (∀ x, 0 ≤ x → 0 < channelE2 x) ∧
  (∀ x, 0 ≤ x → x ≤ (2637 / 100) → 0 < inner2_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ (2637 / 100) → 0 < innerResidual2_0 x) ∧
  (∀ x, (2637 / 100) ≤ x → 0 < tail2 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA3 x) ∧
  (∀ x, 0 ≤ x → 0 < channelE3 x) ∧
  (∀ x, 0 ≤ x → x ≤ (367 / 25) → 0 < inner3_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ (367 / 25) → 0 < innerResidual3_0 x) ∧
  (∀ x, (367 / 25) ≤ x → 0 < tail3 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA4 x) ∧
  (∀ x, 0 ≤ x → 0 < channelE4 x) ∧
  (∀ x, 0 ≤ x → x ≤ (144579 / 25000) → 0 < inner4_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ (144579 / 25000) → 0 < innerResidual4_0 x) ∧
  (∀ x, (144579 / 25000) ≤ x → x ≤ (3047 / 100) → 0 < inner4_1 x) ∧
  (∀ x, (144579 / 25000) ≤ x → x ≤ (3047 / 100) → 0 < innerResidual4_1 x) ∧
  (∀ x, (3047 / 100) ≤ x → 0 < tail4 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA5 x) ∧
  (∀ x, 0 ≤ x → 0 < channelE5 x) ∧
  (∀ x, 0 ≤ x → x ≤ (367 / 25) → 0 < inner5_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ (367 / 25) → 0 < innerResidual5_0 x) ∧
  (∀ x, (367 / 25) ≤ x → 0 < tail5 x) ∧
  (∀ x, 0 ≤ x → 0 < channelA6 x) ∧
  (∀ x, 0 ≤ x → 0 < channelE6 x) ∧
  (∀ x, 0 ≤ x → x ≤ (2637 / 100) → 0 < inner6_0 x) ∧
  (∀ x, 0 ≤ x → x ≤ (2637 / 100) → 0 < innerResidual6_0 x) ∧
  (∀ x, (2637 / 100) ≤ x → 0 < tail6 x) ∧
  (∀ x, 0 ≤ x → x ≤ 1 → 0 < majorant0 x) ∧
  (∀ x, 1 ≤ x → x ≤ 2 → 0 < majorant1 x) ∧
  (∀ x, 2 ≤ x → x ≤ 4 → 0 < majorant2 x) ∧
  (∀ x, 4 ≤ x → x ≤ 8 → 0 < majorant3 x) ∧
  (∀ x, 8 ≤ x → x ≤ 16 → 0 < majorant4 x) ∧
  (∀ x, 16 ≤ x → x ≤ 32 → 0 < majorant5 x) ∧
  (∀ x, 32 ≤ x → x ≤ 64 → 0 < majorant6 x) ∧
  (∀ x, 64 ≤ x → x ≤ 112 → 0 < majorant7 x) ∧
  (∀ x, 112 ≤ x → x ≤ 160 → 0 < majorant8 x)

theorem all_polynomial_conditions : PolynomialConditions :=
  ⟨denominatorU_pos, denominatorL_pos, negativeResidualU_pos, residualL_pos, channelA2_pos, channelE2_pos, inner2_0_pos, innerResidual2_0_pos, tail2_pos, channelA3_pos, channelE3_pos, inner3_0_pos, innerResidual3_0_pos, tail3_pos, channelA4_pos, channelE4_pos, inner4_0_pos, innerResidual4_0_pos, inner4_1_pos, innerResidual4_1_pos, tail4_pos, channelA5_pos, channelE5_pos, inner5_0_pos, innerResidual5_0_pos, tail5_pos, channelA6_pos, channelE6_pos, inner6_0_pos, innerResidual6_0_pos, tail6_pos, majorant0_pos, majorant1_pos, majorant2_pos, majorant3_pos, majorant4_pos, majorant5_pos, majorant6_pos, majorant7_pos, majorant8_pos⟩

end GLVortex
