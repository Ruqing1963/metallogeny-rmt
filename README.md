# metallogeny-rmt

**Charge-and-Release Rhythms in Ore-Forming Systems: A Random Matrix Theory Test of Metallogenic Timing across Four Provinces**

Single ore systems show level repulsion whether or not magmatism is involved.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## Core Finding

The **charge-and-release hypothesis**: each mineralization pulse depletes the
local metal/fluid/heat budget, which must re-accumulate before the next pulse —
imposing level repulsion. Tested across four metallogenic provinces:

| Province | Setting | n | ⟨r⟩ | 95% CI | Poisson excluded? |
|---|---|---|---|---|---|
| Andean porphyry Cu | subduction | 40 | 0.601 | [0.471, 0.670] | ✓ |
| Nanling W-Sn | intracontinental | 36 | 0.574 | [0.488, 0.685] | ✓ |
| Tethyan porphyry Cu | collision | 21 | 0.712 | [0.510, 0.797] | ✓ |
| **Orogenic gold** | **metamorphic (non-magmatic)** | 22 | 0.678 | [0.606, 0.823] | ✓ |

All four single-source configurations show GOE/GUE repulsion; all mixed
configurations collapse toward Poisson (⟨r⟩ = 0.31–0.47).

## The Non-Magmatic Control — Key Result

**Orogenic gold fluids come from metamorphic devolatilization, not magmatism.**
Yet orogenic gold shows repulsion as strong as the magmatic provinces. Because
no magma is involved, this repulsion **cannot be inherited from magmatic-intrusion
timing** — it supports an *intrinsic* charge-release memory in the ore system
itself (repeated fluid focusing along reactivated shear zones).

Grouped by craton, the gold signal is consistent (Yilgarn 0.65, Abitibi 0.77,
West Africa 0.68; Yilgarn-vs-Abitibi KS p=0.83 = no difference).

## Grand Unification — Five Systems

| System | Mechanism | single ⟨r⟩ | mixed | separation |
|---|---|---|---|---|
| Mantle plumes | core-mantle thermal shadow | 0.630 | 0.397 | +0.23 |
| Chile porphyry Cu | subduction magma | 0.601 | 0.470 | +0.13 |
| Nanling W-Sn | intracontinental granite | 0.574 | 0.307 | +0.27 |
| Tethyan porphyry Cu | collisional magma | 0.712 | 0.391 | +0.32 |
| Orogenic gold | metamorphic (non-magmatic) | 0.678 | 0.432 | +0.25 |

**Five physically unrelated systems, one principle: single long-memory sources
repel; superposed sources randomize.**

## Files
- `paper/` — LaTeX + PDF (7 pp.) with embedded figures
- `code/metallogeny_rmt.py` — four-province pooled RMT + bootstrap
- `data/` — 4 CSVs (Chile, Nanling, Tethyan, orogenic gold)
- `figures/` — 2 publication PDFs
- `results/` — JSON output

## Reproduce
```bash
pip install -r requirements.txt
cd code
python metallogeny_rmt.py ../data
```

## Four-Racetrack RMT Program
1. Geological boundaries (Myr) → GOE — [zenodo 20766310](https://zenodo.org/records/20766310)
2. Seismotectonic rhythms → scale-dependent — [zenodo 20768130](https://zenodo.org/records/20768130)
3. Mantle plumes (Gyr) → single-source GOE — [zenodo 20768420](https://zenodo.org/records/20768420)
4. Metallogeny (this work) → single ore system GOE, magmatic & non-magmatic

## Honest Limitations
- Compiled ages are representative literature values, **not a systematic census**
- Sample sizes modest (21–40 pooled intervals/province)
- GOE-vs-GUE not resolved at present sample sizes
- Definitive intrinsic-vs-extrinsic test needs per-district datasets dense in
  **both** ore and magmatic ages — identified as the key next step

## License
MIT. Geochronology compiled from published literature (cited in paper).
