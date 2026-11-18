# Ginzburg–Landau vortex Morse indices, degrees four through eight

Supplement for *Computer-assisted Morse indices of nonmagnetic Ginzburg–Landau vortices of degrees four through eight*, Davis Rens, SuperMatics Inc. Version 1.0.0, Apache-2.0.

The accompanying SuperMatics Technical Report SM-TR-2026-01 is available at [10.5281/zenodo.24326103](https://doi.org/10.5281/zenodo.24326103). It has not been peer reviewed.

Repository: https://github.com/SuperMatics-Inc/gl-vortex-morse-indices

Archival DOI: [10.5281/zenodo.24326105](https://doi.org/10.5281/zenodo.24326105).

## Results

Under the analytic hypotheses of the report, the full real Morse indices of the standard nonmagnetic vortices of degrees 4 through 8 are 12, 20, 30, 36 and 48, respectively. The result also gives bounded-kernel nondegeneracy at degree four. The conjectured upper bound n(n−1) is attained at degrees 4–6 and is strict at degrees 7–8. The calculation stops at degree eight.

Exact rational certificates bound the vortex profiles and give matching upper and lower index counts. The passage from finite matrices to the continuous operator uses the stated spectral and form-domain hypotheses. The coupled variational argument in `degrees5-8/PROOF.md` is a written proof, outside the Lean formalization. The separate FLINT verifier checks the arithmetic with different code, using the same mathematical reductions.

The optional cell checker uses quadrature and finite-element condensation without rigorous error bounds. The figure shows counts from the certificates. There are no measured data.

| Degree | Certified radial counts, ascending m from 2 | Conditional full real index | n(n−1) |
|---|---|---:|---:|
| 4 | 1,1,2,1,1 | 12 | 12 |
| 5 | 1,1,2,2,2,1,1 | 20 | 20 |
| 6 | 1,1,2,2,3,2,2,1,1 | 30 | 30 |
| 7 | 1,1,2,2,2,3,2,2,1,1,1 | 36 | 42 |
| 8 | 1,1,2,2,2,3,3,3,2,2,1,1,1 | 48 | 56 |

Each radial block has two real angular copies. The m=0 and m=1 blocks have no negative directions under the stated hypotheses.

## Reproduce the calculations

Python 3.13 and the exact package versions in `requirements.txt`:

```sh
uv run --python 3.13 --no-project --with-requirements requirements.txt python -B replay_all.py
```

The driver checks `INPUT-HASHES.json` and runs seven subprocesses. It runs the degree-four barrier, upper, lower and sanity checks, followed by the separate FLINT verification for degree four. It then checks all degree-five-to-eight upper and lower bounds with FLINT and runs the negative controls. All seven must pass. Results go to `$SM_OUT` if set, otherwise to `logs/merged-replay/`. GitHub Actions uses the equivalent entry point, `scripts/replay.py`. The checks run without a network inference service.

Other shipped executable scripts, with the same Python dependency prefix:

```sh
python -B scripts/check_remainder.py
python -B degrees5-8/validation/chk_cells.py 7 6 0,2500,5298
python -B figures/plot_summary.py
```

The first command checks the exact cancellation for the exterior O(r^-8) profile remainder. The written proof uses the maximum principle to turn that calculation into a remainder bound. The second command samples three cells. The third regenerates the PDF, a temporary PNG preview and `data/morse-indices.json`; delete the preview after viewing. `data/execution.json` lists the scripts that were run, including the callers of imported helpers. The proof uses the saved rational trial data; the floating-point scripts that designed them are not included.

## Lean scope

Both projects use Lean and Mathlib v4.35.0. For degree four, run `lake update`, `lake exe cache get`, `lake build` and `lake env lean AxiomCheck.lean` from `lean/`. Fetching the cache avoids rebuilding Mathlib. The reported checks used prebuilt dependencies.

| Scope | Machine-checked statements | Not machine-checked |
|---|---|---|
| Degree four, `lean/` | 40 real-polynomial positivity statements, exact Gram signs and central determinant, square-completion and monotonicity algebra, conditional channel bookkeeping | PDE construction, analytic profile and spectral inputs, differential/integral identification, admissibility, trace/variational/Fourier arguments |
| Degrees five to eight, `degrees5-8/lean/` | 15 of 97 files (96 generated modules plus the root import file), checked modules listed below | Remaining modules; real positivity, cell/Gram/pivot construction and linkage, full certificate chain, variational bridge and PDE |

The higher-degree checks cover integer-list arithmetic in the following modules (prefix `GLDegrees.`): `Finite`, `Semantics`, `Inertia`, `Gram`, `Profile5`, `TailProfile5`, `Tail5M2`, `Pivots5M2` through `Pivots5M8`, and `TailProfile6`. To check an individual module, use `lake build GLDegrees.Finite`, for example. The higher-degree root imports the unchecked modules and does not build as supplied.

`GLDegrees.Tail5M2` passed with a 10000 MB kernel limit after exceeding the earlier 5000 MB limit. The higher-degree lakefile uses `-M 10000`. This is a kernel allowance, not a measured peak-memory requirement. The higher-degree invocation below used prebuilt dependencies; the higher-degree project was not checked from a clean clone. From `degrees5-8/lean/`, prepare the dependencies and build `GLDegrees.TailProfile5`, then run:

```sh
lake env lean -M 10000 -j 1 \
  -DmaxHeartbeats=0 -DmaxRecDepth=200000 \
  -o .lake/build/lib/lean/GLDegrees/Tail5M2.olean \
  GLDegrees/Tail5M2.lean
```

## Contents and citation

- `certificates/`: degree-four exact verifiers, imported helpers and saved trial data.
- `degrees5-8/profiles/`, `degrees5-8/certificates/`: exact profile and compressed upper/lower witnesses.
- `degrees5-8/verify_flint.py`, `semantic_mutants.py`, `validation/`: exact and numerical checks.
- `lean/`, `degrees5-8/lean/`: partial formal sources.
- `scripts/`, `data/`, `figures/`: entry points, calculated summaries, execution coverage and figure sources.

The report is distributed separately. Use `CITATION.cff` to cite the software and cite the report separately.

## AI use

AI agents generated substantive proofs, certificates, code, Lean sources and manuscript drafts and edits under Davis Rens's direction and review. No AI agent is an author. Agreement among models is not proof.
