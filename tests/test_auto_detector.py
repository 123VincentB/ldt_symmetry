# -*- coding: utf-8 -*-
"""tests/test_auto_detector.py

Test LdtAutoDetector sur le jeu de 42 fichiers .ldt anonymisés.

Compare le mode détecté automatiquement (LdtAutoDetector) avec le mode
théorique attendu (TestSet.csv → colonne isym_after).

Résultats :
    - Rapport console : OK / KO par fichier
    - Matrice de confusion (isym attendu vs détecté)
    - Diagnostics détaillés pour les cas KO
    - Export CSV : data/output/auto_detector_results.csv

Usage :
    python tests/test_auto_detector.py
    # Depuis la racine du projet (là où se trouve pyproject.toml)
"""

import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Chemins
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_INPUT   = PROJECT_ROOT / "data" / "input"
DATA_OUTPUT  = PROJECT_ROOT / "data" / "output"
TESTSET_CSV  = DATA_INPUT / "TestSet.csv"
RESULTS_CSV  = DATA_OUTPUT / "auto_detector_results.csv"

DATA_OUTPUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
try:
    from pyldt import LdtReader
except ImportError:
    sys.exit("ERREUR : pyldt introuvable. Active le venv et vérifie pip install -e .")

try:
    from ldt_symmetry import LdtAutoDetector
except ImportError:
    sys.exit("ERREUR : ldt_symmetry introuvable. Vérifie pip install -e .")

# ---------------------------------------------------------------------------
# Lecture TestSet.csv
# ---------------------------------------------------------------------------
def load_testset(csv_path: Path) -> list[dict]:
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        return list(reader)


# ---------------------------------------------------------------------------
# Helpers affichage
# ---------------------------------------------------------------------------
ISYM_LABEL = {0: "none", 1: "rot", 2: "C0", 3: "C90", 4: "both"}
COL_W = 12

def _fmt(val, width=COL_W):
    return str(val).ljust(width)

def _scores_str(scores: dict) -> str:
    return "  ".join(f"ISYM{k}={v:.4f}" for k, v in sorted(scores.items()))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    rows = load_testset(TESTSET_CSV)
    detector = LdtAutoDetector()

    results = []       # liste de dict pour export CSV
    errors  = []       # fichiers en erreur (exception)
    ko_list = []       # fichiers KO (détection incorrecte)

    # En-tête tableau console
    header_line = (
        _fmt("Fichier")
        + _fmt("Attendu")
        + _fmt("Détecté")
        + _fmt("Candidat")
        + _fmt("Accepté")
        + "Résultat"
    )
    sep = "-" * len(header_line)
    print()
    print("=" * len(header_line))
    print("  TEST LdtAutoDetector — 42 fichiers")
    print("=" * len(header_line))
    print(header_line)
    print(sep)

    for row in rows:
        name         = row["file"].strip()
        isym_expected = int(row["isym_after"].strip())
        ldt_path      = DATA_INPUT / f"{name}.ldt"

        if not ldt_path.exists():
            print(f"{name}: FICHIER MANQUANT → {ldt_path}")
            errors.append({"file": name, "status": "MISSING"})
            continue

        try:
            ldt = LdtReader.read(ldt_path)
            isym_detected, diag = detector.detect(ldt, return_diag=True)
        except Exception as exc:
            print(f"{_fmt(name)}EXCEPTION : {exc}")
            errors.append({"file": name, "status": f"EXCEPTION: {exc}"})
            continue

        isym_candidate = diag["isym_candidate"]
        accepted       = diag["accepted"]
        scores         = diag["scores"]
        ok             = (isym_detected == isym_expected)
        status         = "OK" if ok else "KO ←"

        line = (
            _fmt(name)
            + _fmt(f"{isym_expected} ({ISYM_LABEL[isym_expected]})")
            + _fmt(f"{isym_detected} ({ISYM_LABEL[isym_detected]})")
            + _fmt(f"{isym_candidate} ({ISYM_LABEL[isym_candidate]})")
            + _fmt(str(accepted))
            + status
        )
        print(line)

        result_row = {
            "file":           name,
            "isym_expected":  isym_expected,
            "isym_detected":  isym_detected,
            "isym_candidate": isym_candidate,
            "accepted":       accepted,
            "score_1":        f"{scores[1]:.6f}",
            "score_2":        f"{scores[2]:.6f}",
            "score_3":        f"{scores[3]:.6f}",
            "score_4":        f"{scores[4]:.6f}",
            "status":         "OK" if ok else "KO",
        }
        results.append(result_row)
        if not ok:
            ko_list.append((name, isym_expected, diag))

    print(sep)

    # ---------------------------------------------------------------------------
    # Résumé global
    # ---------------------------------------------------------------------------
    n_total = len(results)
    n_ok    = sum(1 for r in results if r["status"] == "OK")
    n_ko    = sum(1 for r in results if r["status"] == "KO")
    n_err   = len(errors)

    print()
    print(f"  Total   : {n_total + n_err}")
    print(f"  OK      : {n_ok}")
    print(f"  KO      : {n_ko}")
    print(f"  Erreurs : {n_err}")
    print(f"  Taux    : {100*n_ok/n_total:.1f} %" if n_total else "")
    print()

    # ---------------------------------------------------------------------------
    # Matrice de confusion (attendu = lignes, détecté = colonnes)
    # ---------------------------------------------------------------------------
    isym_labels = [0, 1, 2, 3, 4]
    confusion = {exp: {det: 0 for det in isym_labels} for exp in isym_labels}
    for r in results:
        confusion[r["isym_expected"]][r["isym_detected"]] += 1

    print("  Matrice de confusion (lignes=attendu, colonnes=détecté)")
    print()
    header_conf = "        " + "".join(f"  det={k}".ljust(8) for k in isym_labels)
    print(header_conf)
    print("  " + "-" * (len(header_conf) - 2))
    for exp in isym_labels:
        row_str = f"  exp={exp} |"
        for det in isym_labels:
            val = confusion[exp][det]
            marker = " *" if (exp != det and val > 0) else "  "
            row_str += f"{marker}{val}".ljust(8)
        print(row_str)
    print()
    print("  (* = erreur de classification)")
    print()

    # ---------------------------------------------------------------------------
    # Diagnostics KO
    # ---------------------------------------------------------------------------
    if ko_list:
        print("=" * 60)
        print("  DÉTAIL DES CAS KO")
        print("=" * 60)
        for name, isym_exp, diag in ko_list:
            print()
            print(f"  Fichier   : {name}")
            print(f"  Attendu   : {isym_exp} ({ISYM_LABEL[isym_exp]})")
            print(f"  Détecté   : {diag['isym_detected']} ({ISYM_LABEL[diag['isym_detected']]})")
            print(f"  Candidat  : {diag['isym_candidate']} ({ISYM_LABEL[diag['isym_candidate']]})")
            print(f"  Accepté   : {diag['accepted']}")
            print(f"  Seuils    : score_th={diag['score_threshold']}  rot_score_th={diag['rot_score_threshold']}")
            print(f"  Scores    : {_scores_str(diag['scores'])}")
            shape = diag.get("shape_detail", {})
            if shape:
                print("  Shape detail :")
                for k, v in shape.items():
                    print(f"    {k}: {v}")
            print("  " + "-" * 56)

    # ---------------------------------------------------------------------------
    # Export CSV
    # ---------------------------------------------------------------------------
    if results:
        fieldnames = [
            "file", "isym_expected", "isym_detected", "isym_candidate",
            "accepted", "score_1", "score_2", "score_3", "score_4", "status"
        ]
        with open(RESULTS_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=";")
            writer.writeheader()
            writer.writerows(results)
        print()
        print(f"  Résultats exportés → {RESULTS_CSV}")

    print()
    return n_ko + n_err  # exit code non nul si échecs


if __name__ == "__main__":
    sys.exit(main())
