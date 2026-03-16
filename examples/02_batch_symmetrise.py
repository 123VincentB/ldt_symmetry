# -*- coding: utf-8 -*-
"""examples/02_batch_symmetrise.py

Batch symmetrisation of a folder of EULUMDAT (.ldt) files.

- Detects the ISYM mode automatically for each file.
- Symmetrises files where a plausible mode is found.
- Writes symmetrised files to an output folder.
- Exports a CSV report with ISYM, scores, and status for every file.

Usage
-----
    python examples/02_batch_symmetrise.py

Edit INPUT_DIR, OUTPUT_DIR, and REPORT_CSV below to match your paths.
"""

import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — edit these paths as needed
# ---------------------------------------------------------------------------
INPUT_DIR  = Path("data/input")
OUTPUT_DIR = Path("data/output")
REPORT_CSV = OUTPUT_DIR / "batch_report.csv"

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
try:
    from pyldt import LdtReader, LdtWriter
except ImportError:
    sys.exit("ERROR: pyldt not found. Run: pip install eulumdat-py>=0.1.4")

try:
    from ldt_symmetry import LdtAutoDetector, LdtSymmetriser
except ImportError:
    sys.exit("ERROR: ldt_symmetry not found. Run: pip install eulumdat-symmetry")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
ISYM_LABEL = {0: "none", 1: "rot", 2: "C0", 3: "C90", 4: "both"}


def process_file(
    path: Path,
    detector: LdtAutoDetector,
    output_dir: Path,
) -> dict:
    """Process a single .ldt file. Returns a result dict for the CSV report."""
    result = {
        "file":           path.name,
        "isym_detected":  "",
        "isym_candidate": "",
        "accepted":       "",
        "score_1":        "",
        "score_2":        "",
        "score_3":        "",
        "score_4":        "",
        "output_file":    "",
        "status":         "",
    }

    try:
        ldt = LdtReader.read(path)
    except Exception as exc:
        result["status"] = f"READ ERROR: {exc}"
        return result

    try:
        isym, diag = detector.detect(ldt, return_diag=True)
    except Exception as exc:
        result["status"] = f"DETECT ERROR: {exc}"
        return result

    scores = diag["scores"]
    result.update({
        "isym_detected":  isym,
        "isym_candidate": diag["isym_candidate"],
        "accepted":       diag["accepted"],
        "score_1":        f"{scores[1]:.6f}",
        "score_2":        f"{scores[2]:.6f}",
        "score_3":        f"{scores[3]:.6f}",
        "score_4":        f"{scores[4]:.6f}",
    })

    if isym == 0:
        result["status"] = "SKIPPED (no symmetry detected)"
        return result

    try:
        ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
        out_path = output_dir / (path.stem + f"_ISYM{isym}.ldt")
        LdtWriter.write(ldt_sym, out_path)
        result["output_file"] = out_path.name
        result["status"] = "OK"
    except Exception as exc:
        result["status"] = f"WRITE ERROR: {exc}"

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ldt_files = sorted(INPUT_DIR.glob("*.ldt"))
    if not ldt_files:
        print(f"No .ldt files found in {INPUT_DIR.resolve()}")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    detector = LdtAutoDetector()
    results  = []

    # Console header
    col = 14
    header = (
        f"{'File':<30}"
        f"{'Detected':<{col}}"
        f"{'Candidate':<{col}}"
        f"{'Score C0':<{col}}"
        f"{'Score C90':<{col}}"
        f"Status"
    )
    sep = "-" * len(header)
    print()
    print("=" * len(header))
    print(f"  Batch symmetrisation — {len(ldt_files)} file(s)")
    print("=" * len(header))
    print(header)
    print(sep)

    for path in ldt_files:
        r = process_file(path, detector, OUTPUT_DIR)
        results.append(r)

        isym_det = r["isym_detected"]
        isym_lbl = (
            f"{isym_det} ({ISYM_LABEL[int(isym_det)]})"
            if isinstance(isym_det, int) else "—"
        )
        cand = r["isym_candidate"]
        cand_lbl = (
            f"{cand} ({ISYM_LABEL[int(cand)]})"
            if isinstance(cand, int) else "—"
        )
        s2 = r["score_2"] or "—"
        s3 = r["score_3"] or "—"

        print(
            f"{r['file']:<30}"
            f"{isym_lbl:<{col}}"
            f"{cand_lbl:<{col}}"
            f"{s2:<{col}}"
            f"{s3:<{col}}"
            f"{r['status']}"
        )

    print(sep)

    # Summary
    n_ok      = sum(1 for r in results if r["status"] == "OK")
    n_skipped = sum(1 for r in results if r["status"].startswith("SKIPPED"))
    n_errors  = sum(1 for r in results if "ERROR" in r["status"])

    print()
    print(f"  Total   : {len(results)}")
    print(f"  OK      : {n_ok}")
    print(f"  Skipped : {n_skipped}  (ISYM 0 — no symmetry detected)")
    print(f"  Errors  : {n_errors}")
    print()

    # CSV report
    fieldnames = [
        "file", "isym_detected", "isym_candidate", "accepted",
        "score_1", "score_2", "score_3", "score_4",
        "output_file", "status",
    ]
    with open(REPORT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=";")
        writer.writeheader()
        writer.writerows(results)

    print(f"  Report  : {REPORT_CSV.resolve()}")
    print()

    return n_errors


if __name__ == "__main__":
    sys.exit(main())
