# eulumdat-symmetry

[![PyPI](https://img.shields.io/pypi/v/eulumdat-symmetry)](https://pypi.org/project/eulumdat-symmetry/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/eulumdat-symmetry)](https://pypi.org/project/eulumdat-symmetry/)
[![License: MIT](https://img.shields.io/github/license/123VincentB/eulumdat-symmetry)](https://github.com/123VincentB/eulumdat-symmetry/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19047884.svg)](https://doi.org/10.5281/zenodo.19047884)

An extension to [pyldt](https://pypi.org/project/eulumdat-py/) for symmetrising EULUMDAT (`.ldt`) photometric files.

Developed in an ISO 17025 accredited photometry laboratory.

---

## Overview

EULUMDAT files store luminous intensity distributions as a `[mc × ng]` matrix of C-planes and γ-angles. When a luminaire has a symmetric distribution, the file can declare this via the `ISYM` field (1–4), which allows lighting design software to reconstruct the full distribution from a subset of C-planes.

This library provides two tools:

- **`LdtSymmetriser`** — symmetrise a `pyldt.Ldt` object for a known ISYM mode. Deterministic, no configuration needed.
- **`LdtAutoDetector`** — automatically detect the most plausible ISYM mode from the shape of the distribution. Optional; use only when the mode is not known in advance.

---

## Installation

```bash
pip install eulumdat-symmetry
```

Requires [eulumdat-py](https://pypi.org/project/eulumdat-py/) >= 1.0.0 (`pyldt`).

---

## ISYM codes

| ISYM | Description              | C-planes symmetrised                    |
| ---- | ------------------------ | --------------------------------------- |
| 0    | No symmetry              | —                                       |
| 1    | Full rotational symmetry | All planes averaged to a single profile |
| 2    | Symmetry about C0–C180   | C and 360°−C averaged                   |
| 3    | Symmetry about C90–C270  | C and 180°−C averaged                   |
| 4    | Quadrant symmetry        | All 4 quadrant mirrors averaged         |

---

## Usage

### Known symmetry mode

```python
from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtSymmetriser

ldt = LdtReader.read("luminaire.ldt")
ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2)
LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")
```

`LdtSymmetriser.symmetrise()` returns a new `Ldt` object — the original is never modified.

If the file already declares a non-zero ISYM, the function returns an unchanged copy by default (the distribution has already been symmetrised at measurement time). Use `force=True` to re-symmetrise regardless.

```python
ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2, force=True)
```

---

### Automatic detection

```python
from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

ldt = LdtReader.read("luminaire.ldt")

detector = LdtAutoDetector()
isym = detector.detect(ldt)

if isym != 0:
    ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
    LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")
else:
    print("No plausible symmetry detected — file left unchanged.")
```

---

### Automatic detection with diagnostic

When `return_diag=True`, `detect()` returns a `(isym, diag)` tuple. The `diag` dict contains the shape metrics, symmetry scores, and the veto decision — useful for understanding why a particular mode was accepted or rejected.

```python
isym, diag = detector.detect(ldt, return_diag=True)

print(f"Detected ISYM: {isym}")
print(f"Candidate ISYM (shape stage): {diag['isym_candidate']}")
print(f"Accepted: {diag['accepted']}")
print(f"Scores: {diag['scores']}")
# diag['scores'] = {1: 0.12, 2: 0.08, 3: 0.31, 4: 0.31}
# diag['score_threshold'] = 0.23
```

---

### Batch processing

```python
from pathlib import Path
from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

detector = LdtAutoDetector()

for path in Path("ldt_files").glob("*.ldt"):
    ldt = LdtReader.read(path)
    isym = detector.detect(ldt)
    if isym != 0:
        ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
        out = path.with_stem(path.stem + "_SYM")
        LdtWriter.write(ldt_sym, out)
        print(f"{path.name} → ISYM {isym} → {out.name}")
    else:
        print(f"{path.name} → no symmetry detected")
```

---

## How automatic detection works

Detection runs in two stages:

**Stage 1 — Shape analysis**

The intensity distribution is projected onto the floor (γ < 90°) and ceiling (γ > 90°) planes using the point-source illuminance formula E = I·cos³(γ)/d². Weighted moments (centroid offset, ellipse anisotropy, principal axis orientation) and angular harmonics (H4 for quadrant detection) are computed for each hemisphere. A candidate mode is chosen by combining the decisions from both hemispheres.

**Stage 2 — Score veto**

The relative RMS asymmetry is measured for the candidate mode. If the score exceeds the threshold, the candidate is rejected and ISYM 0 is returned. This prevents the shape stage from proposing a symmetric mode when the distribution is actually irregular.

Both stages must agree for a mode to be accepted.

### Validation

The detection algorithm was validated on a dataset of 42 real manufacturer
LDT files covering all five ISYM modes (0–4), selected for their
representativeness from several hundred laboratory measurements.
No misclassification was observed (42/42 correct).

---

## Adjusting detection thresholds

The default thresholds work well for typical manufacturer files. If you observe incorrect decisions on your dataset, you can adjust them at instantiation:

```python
detector = LdtAutoDetector(
    score_th=0.20,      # stricter veto for ISYM 2/3/4 (default: 0.23)
    rot_score_th=0.15,  # stricter veto for ISYM 1 (default: 0.20)
)
```

See the `LdtAutoDetector` docstring for the full list of parameters.

---

## Project structure

```
eulumdat-symmetry/
├── ldt_symmetry/
│   ├── __init__.py          ← Public API (LdtSymmetriser, LdtAutoDetector)
│   ├── symmetriser.py       ← LdtSymmetriser
│   ├── auto_detector.py     ← LdtAutoDetector
│   └── _geometry.py         ← Internal geometry and scoring functions
```

---

## Relation to pyldt

`eulumdat-symmetry` is an extension to `pyldt`, not a replacement. It operates on `pyldt.Ldt` objects and delegates all file I/O to `pyldt.LdtReader` and `pyldt.LdtWriter`.

| Task                                   | Tool                                    |
| -------------------------------------- | --------------------------------------- |
| Read / write `.ldt` files              | `pyldt` (`LdtReader`, `LdtWriter`)      |
| Edit header fields                     | `pyldt` (`LdtHeader`)                   |
| Symmetrise the intensity distribution  | `eulumdat-symmetry` (`LdtSymmetriser`)  |
| Detect the symmetry mode automatically | `eulumdat-symmetry` (`LdtAutoDetector`) |

---

## eulumdat-* ecosystem

> **New to the ecosystem?** [eulumdat-quickstart](https://github.com/123VincentB/eulumdat-quickstart) — a step-by-step guide covering all 8 packages with working examples.

| Package | Description |
|---|---|
| [eulumdat-py](https://pypi.org/project/eulumdat-py/) | Read / write EULUMDAT files |
| **`eulumdat-symmetry`** | **Symmetrise and detect ISYM — this package** |
| [eulumdat-plot](https://pypi.org/project/eulumdat-plot/) | Polar intensity diagram (SVG/PNG) |
| [eulumdat-luminance](https://pypi.org/project/eulumdat-luminance/) | Luminance table and polar diagram |
| [eulumdat-ugr](https://pypi.org/project/eulumdat-ugr/) | UGR catalogue (CIE 117/190) |
| [eulumdat-analysis](https://pypi.org/project/eulumdat-analysis/) | Beam half-angle, FWHM |
| [eulumdat-report](https://pypi.org/project/eulumdat-report/) | Full photometric datasheet (HTML/PDF) |
| [eulumdat-ies](https://pypi.org/project/eulumdat-ies/) | LDT ↔ IES LM-63-2002 conversion |

---

## Source code

[github.com/123VincentB/eulumdat-symmetry](https://github.com/123VincentB/eulumdat-symmetry)

---

## License

MIT
