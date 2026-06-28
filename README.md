# metallogeny-rmt

**Charge-and-Release Rhythms in Ore-Forming Systems Are Scalar Clocks, Not Level Repulsion: A Cautionary Random-Matrix Analysis of Metallogenic Timing**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)
[![Version 2.0.0](https://img.shields.io/badge/version-2.0.0-orange.svg)](CHANGELOG.md)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

> ## ⚠️ Version 2.0 — Erratum / self-correction
>
> **This release retracts the central claim of v1.** Version 1 reported that
> single ore-forming systems show **GOE/GUE level repulsion** in pulse timing and
> read it as an intrinsic "charge-and-release memory." **That interpretation was
> wrong and is withdrawn.** A high spacing ratio ⟨r⟩ is *two-ended*: it is equally
> produced by a narrow-marginal **scalar renewal clock** that carries no
> sequential repulsion at all. **Retained as real:** superposition of independent
> sources drives ⟨r⟩ toward the Poisson value. Version 1 is preserved (Git tag
> `v1.0`, `paper/archive/v1/`, and as a prior version of the Zenodo record). See
> [`CHANGELOG.md`](CHANGELOG.md).

## What changed and why

A high nearest-neighbour spacing ratio (⟨r⟩ ≈ 0.5–0.7) does **not** by itself
imply level repulsion. It has two distinct sources:

1. **Sequential level repulsion** (the RMT/GOE reading v1 assumed) — neighbouring
   events *avoid* each other in time.
2. **A narrow interval marginal** — a **scalar renewal clock** (relaxation
   oscillator) whose intervals simply cluster around a mean, with *no* ordering
   structure, produces the *same* high ⟨r⟩.

Four diagnostics separate these, and all point to (2) for single ore systems:

| Diagnostic | Result | Reading |
|---|---|---|
| **Shuffle test** (reorder intervals) | ⟨r⟩ barely moves (gold obs−shuffle = −0.031) | signal is in the *marginal*, not the *sequence* |
| **Zero-repulsion gamma clock** (k = 1/CV², finite-N matched) | brackets every observed ⟨r⟩ (gold 0.726 [0.633, 0.808] ⊃ 0.678) | a no-repulsion clock reproduces the data |
| **Bohigas–Giannoni–Schmit** | a scalar charge-release reservoir is integrable | **structurally not GOE** |
| **Closed-data control** (orogenic gold) | zero-repulsion null through v1 normalisation gives obs−shuffle = −0.028 ± 0.021 | apparent rigidity is a compositional artifact |

For the non-magmatic orogenic-gold "control," the within-deposit intervals are
differences between *independent geochronometers on the same deposit* (e.g. U-Pb,
Ar-Ar, Re-Os), so the spacings are dominated by analytical age errors and the
v1 within-deposit normalisation (Σsᵢ = const) forces the apparent rigidity.

## What is retained (still real)

**Superposition → Poisson.** Pooling independent ore sources drives the spacing
ratio toward the Poisson value (⟨r⟩ → 0.386). This survives the audit unchanged.

| Configuration | ⟨r⟩ |
|---|---|
| Single source (narrow scalar clock) | high (0.57–0.71), **no repulsion** |
| Mixed / pooled independent sources | → 0.386 (Poisson) |

## Reference values

Poisson ⟨r⟩ = 0.386 · GOE ⟨r⟩ = 0.531 · GUE ⟨r⟩ = 0.603

## Files
- `paper/` — LaTeX + PDF (v2.0 erratum) with embedded audit figures
- `paper/archive/v1/` — preserved v1 paper and figures
- `code/audit_scalar_clock.py` — four-panel audit (shuffle, gamma clock, closed-data control, age-error MC)
- `code/make_figs_v2.py` — regenerates the v2 figures
- `code/metallogeny_rmt.py` — original four-province pooled analysis (now also reports shuffle test + CV)
- `data/` — 4 CSVs (Chile, Nanling, Tethyan, orogenic gold)
- `figures/` — publication PDFs
- `results/` — JSON output (`four_province_results.json`, `audit_results.json`)

## Reproduce
```bash
pip install -r requirements.txt
cd code
python metallogeny_rmt.py ../data      # reproduces v1 ⟨r⟩ + shuffle test + CV
python audit_scalar_clock.py           # four-panel scalar-clock audit
python make_figs_v2.py                 # regenerates v2 figures
```

## Method basis
Scalar-clock audit method: Zenodo DOI [10.5281/zenodo.20980670](https://doi.org/10.5281/zenodo.20980670).

## On self-correction
This erratum is issued as a Zenodo **new version** with v1 fully preserved and
citable. Transparent correction — not silent revision — is the intended standard
for this work.

## License
MIT. Geochronology compiled from published literature (cited in paper).
