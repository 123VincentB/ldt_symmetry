# -*- coding: utf-8 -*-
"""01_basic_usage.py
=====================
Basic usage examples for eulumdat-symmetry.

Demonstrates:
  1. Symmetrise with a known ISYM mode
  2. Automatic detection of the symmetry mode
  3. Automatic detection with diagnostic information

Requirements:
    pip install eulumdat-symmetry
"""

from pathlib import Path

from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

# ---------------------------------------------------------------------------
# Paths — adjust to your own .ldt file
# ---------------------------------------------------------------------------
INPUT_FILE = Path("luminaire.ldt")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Example 1 — Symmetrise with a known ISYM mode
# ---------------------------------------------------------------------------
# Use this when you already know the symmetry of the luminaire.
# ISYM codes:
#   0 = no symmetry
#   1 = full rotational symmetry
#   2 = symmetry about C0-C180
#   3 = symmetry about C90-C270
#   4 = quadrant symmetry (C0-C180 and C90-C270)

ldt = LdtReader.read(INPUT_FILE)
ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2)
LdtWriter.write(ldt_sym, OUTPUT_DIR / "luminaire_ISYM2.ldt")

print(f"Example 1 — ISYM after symmetrisation: {ldt_sym.header.isym}")

# ---------------------------------------------------------------------------
# Example 2 — Automatic detection
# ---------------------------------------------------------------------------
# Use this when the symmetry mode is not known in advance.
# Returns 0 if no plausible symmetry is detected.

ldt = LdtReader.read(INPUT_FILE)
detector = LdtAutoDetector()
isym = detector.detect(ldt)

if isym != 0:
    ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
    LdtWriter.write(ldt_sym, OUTPUT_DIR / f"luminaire_ISYM{isym}.ldt")
    print(f"Example 2 — Detected ISYM: {isym}")
else:
    print("Example 2 — No plausible symmetry detected — file left unchanged.")

# ---------------------------------------------------------------------------
# Example 3 — Automatic detection with diagnostic information
# ---------------------------------------------------------------------------
# return_diag=True returns a (isym, diag) tuple with internal metrics.
# Useful for understanding why a particular mode was accepted or rejected.

ldt = LdtReader.read(INPUT_FILE)
isym, diag = detector.detect(ldt, return_diag=True)

print(f"\nExample 3 — Diagnostic:")
print(f"  Detected ISYM  : {diag['isym_detected']}")
print(f"  Candidate ISYM : {diag['isym_candidate']}  (shape stage)")
print(f"  Accepted       : {diag['accepted']}")
print(f"  Scores         : {diag['scores']}")
print(f"  Threshold      : {diag['score_threshold']}  (ISYM 2/3/4)")
print(f"  Threshold rot  : {diag['rot_score_threshold']}  (ISYM 1)")
