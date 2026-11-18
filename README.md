# Ginzburg–Landau vortex Morse indices, degrees four through eight

Supplement for *Computer-assisted Morse indices of nonmagnetic Ginzburg–Landau vortices of degrees four through eight*, Davis Rens, SuperMatics Inc. Version 1.0.0, Apache-2.0.

Accompanying report: **SM-TR-2026-01**, technical report, not peer reviewed. Report DOI: [10.5281/zenodo.24326103](https://doi.org/10.5281/zenodo.24326103).

Repository: https://github.com/SuperMatics-Inc/gl-vortex-morse-indices

Archival DOI: [10.5281/zenodo.24326105](https://doi.org/10.5281/zenodo.24326105).

## Results and boundaries

**PROVED, conditional on the analytic hypotheses stated in the manuscript:** the full real Morse indices of the standard nonmagnetic vortices of degrees 4, 5, 6, 7 and 8 are respectively 12, 20, 30, 36 and 48. Bounded-kernel nondegeneracy is proved for degree four only. The conjectured upper bound n(n−1) is attained at degrees 4–6 and strictly non-sharp at degrees 7–8; this does not refute that upper bound. No assertion is made for degree nine or above.

**CERTIFIED:** exact-rational finite profile, negative-subspace and upper-index certificates. Their passage to the continuous differential operator uses the written analytic arguments and imported spectral and form-domain hypotheses. `degrees5-8/PROOF.md` supplies the coupled variational argument; it is not formalised in Lean. Code independence of the FLINT implementation is not independent mathematical derivation. Model agreement is not proof.

**NUMERICAL:** the optional cell checker uses quadrature and finite-element condensation without rigorous error bounds. It is not a certificate. The summary figure displays certificate counts, not new numerical spectra. There are no measured data.

| Degree | Certified radial counts, ascending m from 2 | Conditional full real index | n(n−1) |
|---|---|---:|---:|
| 4 | 1,1,2,1,1 | 12 | 12 |
| 5 | 1,1,2,2,2,1,1 | 20 | 20 |
| 6 | 1,1,2,2,3,2,2,1,1 | 30 | 30 |
| 7 | 1,1,2,2,2,3,2,2,1,1,1 | 36 | 42 |
| 8 | 1,1,2,2,2,3,3,3,2,2,1,1,1 | 48 | 56 |

Each radial block has two real angular copies. The m=0 and m=1 blocks have no negative directions under the stated hypotheses.

## Replay

Python 3.13 and the exact package versions in `requirements.txt`:

```sh
uv run --python 3.13 --no-project --with-requirements requirements.txt python -B replay_all.py
```

The runner checks `INPUT-HASHES.json` and executes seven children: the degree-four barrier, upper, lower and sanity checks; the code-independent degree-four FLINT check; all degree-five-to-eight FLINT sandwiches; and semantic adversaries. Any failing child makes the runner fail. Results go to `$SM_OUT` if set, otherwise `logs/merged-replay/`. `scripts/replay.py` is the equivalent CI entry point. No network inference service is used.

Other shipped executable scripts, with the same Python dependency prefix:

```sh
python -B scripts/check_remainder.py
python -B degrees5-8/validation/chk_cells.py 7 6 0,2500,5298
python -B figures/plot_summary.py
```

The first checks the exact cancellation algebra for the exterior O(r^-8) profile remainder. The maximum-principle proof is mathematical, not established by symbolic cancellation alone. The second samples three cells only, not every cell or degree. The third regenerates the PDF, a temporary PNG preview and `data/morse-indices.json`; delete the preview after viewing. `data/execution.json` records exercised script coverage. Imported helpers are exercised by their listed callers. Floating-point witness designers are not distributed; the proof uses the saved rationalised data.

## Lean scope

Both projects pin Lean/Mathlib v4.35.0. Obtain Mathlib's prebuilt cache rather than rebuilding Mathlib. For degree four, from `lean/`, the intended commands are `lake update`, `lake exe cache get`, `lake build`, and `lake env lean AxiomCheck.lean`. The portable Git-based dependency configuration and CI have not been tested in a clean clone; the checked results used prebuilt dependencies.

| Scope | Machine-checked statements | Not machine-checked |
|---|---|---|
| Degree four, `lean/` | 40 real-polynomial positivity statements, exact Gram signs and central determinant, square-completion and monotonicity algebra, conditional channel bookkeeping | PDE construction, analytic profile and spectral inputs, differential/integral identification, admissibility, trace/variational/Fourier arguments |
| Degrees five to eight, `degrees5-8/lean/` | 15 of 97 files (96 generated modules plus the root import file), checked modules listed below | Remaining modules; real positivity, cell/Gram/pivot construction and linkage, full certificate chain, variational bridge and PDE |

Checked higher-degree modules (prefix `GLDegrees.`): `Finite`, `Semantics`, `Inertia`, `Gram`, `Profile5`, `TailProfile5`, `Tail5M2`, `Pivots5M2` through `Pivots5M8`, and `TailProfile6`. To check an individual module after building its dependencies, use `lake build GLDegrees.Finite`, for example. The higher-degree root imports unchecked modules and does not build as shipped. No complete degree-five or higher Lean certificate-chain theorem is claimed. The 15-module statement is limited to integer-list arithmetic, not the analytic theorem.

The denominator 97 counts 96 generated modules plus the root import file. `GLDegrees.Tail5M2` was accepted only with a 10000 MB kernel limit; it exceeded the earlier 5000 MB limit. The shipped higher-degree lakefile now sets `-M 10000` to match that accepted check. This is a kernel allowance, not a measured peak-memory requirement. From `degrees5-8/lean/`, after preparing the dependencies and building `GLDegrees.TailProfile5`, the exact command corresponding to that check is:

```sh
lake env lean -M 10000 -j 1 \
  -DmaxHeartbeats=0 -DmaxRecDepth=200000 \
  -o .lake/build/lib/lean/GLDegrees/Tail5M2.olean \
  GLDegrees/Tail5M2.lean
```

The explicit kernel allowance matches the shipped higher-degree setting. This command reproduces the recorded invocation, not an additional clean-clone build: the portable Git-based Lean configuration remains untested.

## Layout and attribution

- `certificates/`: degree-four exact verifiers, imported helpers and saved trial data.
- `degrees5-8/profiles/`, `degrees5-8/certificates/`: exact profile and compressed upper/lower witnesses.
- `degrees5-8/verify_flint.py`, `semantic_mutants.py`, `validation/`: exact and numerical checks.
- `lean/`, `degrees5-8/lean/`: partial formal sources.
- `scripts/`, `data/`, `figures/`: entry points, calculated summaries, execution coverage and figure sources.

The manuscript is distributed separately, not in this software tree. AI agents generated substantive proofs, certificates, implementations, Lean sources and manuscript drafts and edits under the author's direction and review. No LLM is an author. Cite the software using `CITATION.cff` and cite the accompanying manuscript separately.
