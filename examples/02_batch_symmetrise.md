# Batch symmetrisation

This example shows how to process a folder of EULUMDAT files in one pass:
automatic detection, symmetrisation, and a CSV report of results.

## Prerequisites

```bash
pip install eulumdat-symmetry
```

Requires [eulumdat-py](https://pypi.org/project/eulumdat-py/) >= 1.0.0 (`pyldt`).

---

## Script

[`02_batch_symmetrise.py`](02_batch_symmetrise.py)

---

## Configuration

Edit the three path constants at the top of the script:

```python
INPUT_DIR  = Path("data/input")   # folder containing .ldt files
OUTPUT_DIR = Path("data/output")  # symmetrised files written here
REPORT_CSV = OUTPUT_DIR / "batch_report.csv"
```

Then run from the project root:

```bash
python examples/02_batch_symmetrise.py
```

---

## What the script does

For each `.ldt` file in `INPUT_DIR`:

1. **Reads** the file with `LdtReader`.
2. **Detects** the ISYM mode automatically with `LdtAutoDetector`.
3. **Symmetrises** the distribution with `LdtSymmetriser` if a plausible
   mode is found (ISYM 1–4).  Files returning ISYM 0 are skipped — no
   plausible symmetry was detected for that distribution.
4. **Writes** the symmetrised file to `OUTPUT_DIR` with the suffix
   `_ISYM{n}` added to the original stem.
5. **Reports** every file to the console and to `batch_report.csv`.

---

## Console output

```
================================================================
  Batch symmetrisation — 42 file(s)
================================================================
File                          Detected      Candidate     Score C0      Score C90     Status
------------------------------------------------------------------------------------------------
sample_01.ldt                 3 (C90)       3 (C90)       0.412748      0.001823      OK
sample_02.ldt                 4 (both)      4 (both)      0.003201      0.004870      OK
...
sample_24.ldt                 0 (none)      3 (C90)       0.412748      0.301562      SKIPPED (no symmetry detected)
...
------------------------------------------------------------------------------------------------
  Total   : 42
  OK      : 38
  Skipped : 4  (ISYM 0 — no symmetry detected)
  Errors  : 0

  Report  : C:\...\data\output\batch_report.csv
```

**Score columns** show the relative RMS asymmetry for the two most common
mirror planes (C0–C180 and C90–C270).  A low score (< 0.23) means the
distribution is consistent with that symmetry.  See the
[diagnostic reference](01_basic_usage.md#reading-the-diagnostic) for details.

---

## CSV report

`batch_report.csv` contains one row per file:

| Column           | Description                                              |
| ---------------- | -------------------------------------------------------- |
| `file`           | Input filename                                           |
| `isym_detected`  | Final detected ISYM (0–4)                                |
| `isym_candidate` | ISYM proposed by the shape stage (before score veto)     |
| `accepted`       | `True` if the candidate passed the score veto            |
| `score_1`        | Relative RMS score for ISYM 1 (rotational)               |
| `score_2`        | Relative RMS score for ISYM 2 (C0–C180)                  |
| `score_3`        | Relative RMS score for ISYM 3 (C90–C270)                 |
| `score_4`        | Relative RMS score for ISYM 4 (quadrant)                 |
| `output_file`    | Name of the written output file (empty if skipped)       |
| `status`         | `OK`, `SKIPPED (no symmetry detected)`, or error message |

---

## Output filenames

Symmetrised files are written with `_ISYM{n}` appended to the original stem:

```
sample_01.ldt   →   sample_01_ISYM3.ldt
sample_02.ldt   →   sample_02_ISYM4.ldt
```

---

## Error handling

Read errors, detection errors, and write errors are caught per file and
reported in the `status` column of the CSV.  A non-zero exit code is
returned if any errors occurred, which makes the script compatible with
CI pipelines.

---

## See also

- [Basic usage](01_basic_usage.md) — single-file examples and diagnostic reference
- [pyldt documentation](https://pypi.org/project/eulumdat-py/) — read/write EULUMDAT files
- [eulumdat-symmetry on PyPI](https://pypi.org/project/eulumdat-symmetry/)
