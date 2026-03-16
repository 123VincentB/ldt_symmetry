# Basic usage

This example covers the three main use cases of `eulumdat-symmetry`.

## Prerequisites

```bash
pip install eulumdat-symmetry
```

Requires [eulumdat-py](https://pypi.org/project/eulumdat-py/) >= 0.1.4 (`pyldt`).

---

## Example 1 — Symmetrise with a known ISYM mode

Use this when you already know the symmetry of the luminaire.

```python
from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtSymmetriser

ldt = LdtReader.read("luminaire.ldt")
ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2)
LdtWriter.write(ldt_sym, "luminaire_ISYM2.ldt")
```

`symmetrise()` returns a new `Ldt` object — the original is never modified.

### ISYM codes

| ISYM | Symmetry                                 |
| ---- | ---------------------------------------- |
| 0    | No symmetry (identity)                   |
| 1    | Full rotational symmetry                 |
| 2    | Symmetry about C0–C180                   |
| 3    | Symmetry about C90–C270                  |
| 4    | Quadrant symmetry (C0–C180 and C90–C270) |

### Re-symmetrising an already symmetrised file

If the file already declares a non-zero ISYM, `symmetrise()` returns an
unchanged copy by default. Use `force=True` to re-symmetrise regardless:

```python
ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2, force=True)
```

---

## Example 2 — Automatic detection

Use this when the symmetry mode is not known in advance.

```python
from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

ldt = LdtReader.read("luminaire.ldt")

detector = LdtAutoDetector()
isym = detector.detect(ldt)

if isym != 0:
    ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
    LdtWriter.write(ldt_sym, f"luminaire_ISYM{isym}.ldt")
else:
    print("No plausible symmetry detected — file left unchanged.")
```

Detection returns `0` if no symmetric mode is plausible for the distribution.

---

## Example 3 — Automatic detection with diagnostic information

`return_diag=True` returns a `(isym, diag)` tuple with internal metrics.
Useful for understanding why a particular mode was accepted or rejected.

```python
isym, diag = detector.detect(ldt, return_diag=True)

print(f"Detected ISYM  : {diag['isym_detected']}")
print(f"Candidate ISYM : {diag['isym_candidate']}  (shape stage)")
print(f"Accepted       : {diag['accepted']}")
print(f"Scores         : {diag['scores']}")
print(f"Threshold      : {diag['score_threshold']}  (ISYM 2/3/4)")
print(f"Threshold rot  : {diag['rot_score_threshold']}  (ISYM 1)")
```

### Reading the diagnostic

| Key                   | Description                                                                 |
| --------------------- | --------------------------------------------------------------------------- |
| `isym_detected`       | Final result (0–4)                                                          |
| `isym_candidate`      | Mode proposed by the shape analysis stage                                   |
| `accepted`            | `True` if the candidate passed the score veto                               |
| `scores`              | Relative RMS asymmetry per mode `{1: float, 2: float, 3: float, 4: float}` |
| `score_threshold`     | Veto threshold for ISYM 2/3/4 (default: 0.23)                              |
| `rot_score_threshold` | Veto threshold for ISYM 1 (default: 0.20)                                  |

A score close to 0 means the distribution is highly symmetric for that mode.
A score above the threshold means the candidate was rejected — ISYM 0 is returned.

---

## See also

- [pyldt documentation](https://pypi.org/project/eulumdat-py/) — read/write EULUMDAT files
- [eulumdat-symmetry on PyPI](https://pypi.org/project/eulumdat-symmetry/)
- [Source code](https://github.com/123VincentB/ldt_symmetry)
