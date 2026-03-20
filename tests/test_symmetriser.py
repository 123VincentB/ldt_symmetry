# -*- coding: utf-8 -*-
"""test_symmetriser.py
======================
Test de LdtSymmetriser sur le jeu de fichiers réels défini dans TestSet.csv.

Pour chaque fichier :
  - lecture via LdtReader
  - symétrisation avec l'ISYM attendu (depuis le CSV)
  - vérification que header.isym est bien mis à jour
  - écriture du fichier symétrisé dans data/output/

Usage :
    python tests/test_symmetriser.py
    pytest tests/test_symmetriser.py
"""

import csv
from pathlib import Path

from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtSymmetriser

# --- Chemins ---------------------------------------------------------------
ROOT = Path(__file__).parent.parent
INPUT_DIR = ROOT / "data/input"
OUTPUT_DIR = ROOT / "data/output"
CSV_PATH = ROOT / "data/input/TestSet.csv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Test ------------------------------------------------------------------

def test_symmetriser():
    print(f"{'Fichier':<45} {'ISYM':>5}  {'Résultat'}")
    print("-" * 75)

    ok = 0
    skipped = 0
    errors = 0

    with open(CSV_PATH, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            base = row["file"]
            isym_target = int(row["isym_after"])
            ldt_path = INPUT_DIR / f"{base}.ldt"

            if not ldt_path.exists():
                print(f"{'[SKIP] ' + base:<45} {'':>5}  fichier introuvable")
                skipped += 1
                continue

            try:
                ldt = LdtReader.read(ldt_path)
                ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym_target, force=True)

                isym_result = int(ldt_sym.header.isym)
                status = "OK" if isym_result == isym_target else f"ERREUR (isym={isym_result})"

                out_path = OUTPUT_DIR / f"{base}_SYM{isym_target}.ldt"
                LdtWriter.write(ldt_sym, out_path)

                print(f"{base:<45} {isym_target:>5}  {status}")
                ok += 1

            except Exception as ex:
                print(f"{'[ERR] ' + base:<45} {isym_target:>5}  {ex}")
                errors += 1

    print("-" * 75)
    print(f"OK: {ok}  |  Skipped: {skipped}  |  Erreurs: {errors}")

    assert errors == 0, f"{errors} fichier(s) en erreur"


if __name__ == "__main__":
    test_symmetriser()
