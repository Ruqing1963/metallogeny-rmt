# Changelog

All notable changes to this project are documented here. This project preserves
every released version; corrections are issued transparently as new versions
rather than by silent revision.

## [2.0.0] - 2026-06-27 — Erratum / self-correction

This is a **self-correcting erratum**. It retracts the central claim of v1 and
reframes the phenomenon. Version 1 is preserved (Git tag `v1.0`, archived under
`paper/archive/v1/`, and as a prior version of the Zenodo record).

### Retracted
- The claim that **single ore-forming systems exhibit GOE/GUE-class level
  repulsion** in pulse timing, and its interpretation as an intrinsic
  "charge-and-release memory."
- The cross-domain "grand unification of five deep-Earth systems" framing built
  on that claim.
- The use of non-magmatic orogenic gold as a "control" proving intrinsic
  repulsion.

### Why
- A high nearest-neighbour spacing ratio (⟨r⟩ ≈ 0.5–0.7) is **two-ended**: it is
  also produced by a narrow-marginal **scalar renewal clock** (relaxation
  oscillator) that carries **no sequential level repulsion**.
- **Shuffle test:** randomly reordering each province's intervals barely changes
  ⟨r⟩ (gold obs−shuffle = −0.031), so the signal lives in the interval
  *marginal*, not the *sequence*.
- **Zero-repulsion null:** an i.i.d. gamma clock with shape k = 1/CV², matched to
  finite N, **brackets every observed ⟨r⟩** (e.g. gold gamma median 0.726,
  95% [0.633, 0.808], contains observed 0.678).
- **Bohigas–Giannoni–Schmit:** a scalar charge-release reservoir is an integrable
  one-degree-of-freedom system and is therefore **structurally not GOE**.
- **Closed-data artifact:** for orogenic gold the within-deposit "intervals" are
  differences between independent geochronometers on the same deposit; pushing a
  zero-repulsion null through v1's own within-deposit normalisation reproduces the
  apparent rigidity (obs−shuffle = −0.028 ± 0.021), confirming a compositional
  (Σsᵢ = const) artifact.
- **Age-error Monte Carlo:** jittering ages within quoted 1σ spans ⟨r⟩ from
  Poisson to GUE; v1's reported gold value sits at the top edge of that range.

### Retained
- **Superposition → Poisson.** Pooling independent ore sources drives the spacing
  ratio toward the Poisson value (⟨r⟩ → 0.386). This result is real and is kept.

### Added
- `code/audit_scalar_clock.py` — four-panel audit (shuffle test, zero-repulsion
  gamma clock, closed-data control, age-error Monte Carlo).
- `code/make_figs_v2.py` — regenerates figures, with a finite-N 95% band added to
  the master curve for honesty.
- `paper/figs/fig1_scalar_clock.pdf`, `paper/figs/fig2_age_error_superposition.pdf`.
- `results/audit_results.json`.

### Changed
- `paper/paper.tex` fully rewritten as a v2.0 erratum (new title, "Changes from
  version 1" section, reframed abstract, BGS argument, audit methods/results,
  retained superposition section).
- `code/metallogeny_rmt.py` now reports the shuffle test and coefficient of
  variation alongside ⟨r⟩; the misleading "intrinsic memory" conclusion line was
  replaced with the scalar-clock conclusion. v1 ⟨r⟩ values reproduce exactly.
- `README.md` rewritten with a v2 banner and the corrected interpretation.

### Preserved
- v1 paper and figures under `paper/archive/v1/`.
- v1 remains citable as a prior version of the Zenodo record.

### Basis
- Scalar-clock audit method: Zenodo DOI [10.5281/zenodo.20980670](https://doi.org/10.5281/zenodo.20980670).

## [1.0.0] - 2026-06-19

- Initial release. Random-matrix analysis of ore-deposit timing across four
  metallogenic provinces (Andean, Nanling, Tethyan, orogenic gold). **The central
  level-repulsion claim of this version is retracted in v2.0.0; see above.**
