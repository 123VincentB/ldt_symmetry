# -*- coding: utf-8 -*-
"""test_re_symmetrising.py
===========================
Vérifie que la symétrisation est idempotente :
appliquer deux fois le même mode doit donner des valeurs identiques.

Principe :
  1. Lire un fichier ISYM=0 (valeurs brutes)
  2. Symétriser → ldt_sym1
  3. Écrire ldt_sym1, relire → ldt_reread
  4. Re-symétriser avec force=True → ldt_sym2
  5. Comparer ldt_sym1.intensities vs ldt_sym2.intensities

Les valeurs doivent être identiques à la précision numérique près.

Requiert : eulumdat-py >= 0.1.4

Usage :
    python tests/test_re_symmetrising.py
"""

from pathlib import Path

from pyldt import LdtReader, LdtWriter
from ldt_symmetry import LdtSymmetriser

# --- Chemins ---------------------------------------------------------------
ROOT = Path(__file__).parent.parent
INPUT_DIR = ROOT / "data" / "input"
OUTPUT_DIR = ROOT / "data" / "output"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Helpers ---------------------------------------------------------------

def max_abs_diff(mat1, mat2) -> float:
    """Différence absolue maximale entre deux matrices [mc × ng]."""
    return max(
        abs(mat1[i][j] - mat2[i][j])
        for i in range(len(mat1))
        for j in range(len(mat1[i]))
    )

# --- Test ------------------------------------------------------------------

isym0_files = [
    p for p in sorted(INPUT_DIR.glob("*.ldt"))
    if int(LdtReader.read(p).header.isym) == 0
]

if not isym0_files:
    print("Aucun fichier ISYM=0 trouvé dans data/input/ — test ignoré.")
    exit(0)

print(f"{len(isym0_files)} fichier(s) ISYM=0 trouvé(s).\n")
print(f"{'Fichier':<35} {'ISYM':>5}  {'Max diff':>12}  {'Résultat'}")
print("-" * 70)

ok = 0
errors = 0

for ldt_path in isym0_files:
    ldt_orig = LdtReader.read(ldt_path)

    for isym_target in [1, 2, 3, 4]:
        try:
            # Passe 1 : symétriser
            ldt_sym1 = LdtSymmetriser.symmetrise(ldt_orig, isym=isym_target, force=True)

            # Écriture / relecture
            tmp_path = OUTPUT_DIR / f"_tmp_{ldt_path.stem}_isym{isym_target}.ldt"
            LdtWriter.write(ldt_sym1, tmp_path)
            ldt_reread = LdtReader.read(tmp_path)
            tmp_path.unlink()

            # Passe 2 : re-symétriser
            ldt_sym2 = LdtSymmetriser.symmetrise(ldt_reread, isym=isym_target, force=True)

            # Comparaison
            diff = max_abs_diff(ldt_sym1.intensities, ldt_sym2.intensities)
            identical = diff <= 1e-6
            status = "OK" if identical else "DIFF > tolérance"

            print(f"{ldt_path.stem:<35} {isym_target:>5}  {diff:>12.2e}  {status}")

            if identical:
                ok += 1
            else:
                errors += 1

        except Exception as ex:
            print(f"{ldt_path.stem:<35} {isym_target:>5}  {'':>12}  ERREUR: {ex}")
            errors += 1

# --- Résumé ----------------------------------------------------------------
print("-" * 70)
print(f"OK: {ok}  |  Erreurs: {errors}")
