# -*- coding: utf-8 -*-
"""ldt_symmetry.symmetriser
===========================

Provides :class:`LdtSymmetriser`, the main public class for symmetrising
EULUMDAT intensity distributions.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any, List

from ._geometry import safe_angles_deg

# ISYM → internal mode string (used only inside this module)
_ISYM_TO_MODE = {0: "none", 1: "rot", 2: "C0", 3: "C90", 4: "both"}
_VALID_ISYM = frozenset(_ISYM_TO_MODE)


class LdtSymmetriser:
    """Symmetrise a pyldt ``Ldt`` object according to an ISYM mode.

    All methods are static — no instance state is needed.

    ISYM codes
    ----------
    0 : no symmetry (identity — returns a copy with no averaging)
    1 : full rotational symmetry
    2 : symmetry about C0–C180 plane
    3 : symmetry about C90–C270 plane
    4 : quadrant symmetry (C0–C180 and C90–C270)

    Usage
    -----
    ::

        from pyldt import LdtReader, LdtWriter
        from ldt_symmetry import LdtSymmetriser

        ldt = LdtReader.read("luminaire.ldt")
        ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=2)
        LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")
    """

    @staticmethod
    def symmetrise(ldt: Any, isym: int, *, force: bool = False) -> Any:
        """Return a new ``Ldt`` object with the intensity distribution symmetrised.

        The input ``ldt`` is never modified.

        Parameters
        ----------
        ldt:
            A ``pyldt.Ldt`` object as returned by ``LdtReader.read()``.
        isym:
            Target symmetry mode (0–4).  See class docstring.
        force:
            If ``False`` (default) and the file already declares a non-zero
            ISYM, the function returns an unchanged copy — the distribution
            has already been symmetrised at measurement time.
            Set ``force=True`` to re-symmetrise regardless.

        Returns
        -------
        Ldt
            New ``Ldt`` object with updated ``header.isym`` and symmetrised
            ``intensities``.  Header fields other than ``isym`` are unchanged.

        Raises
        ------
        ValueError
            If ``isym`` is not in 0–4.
        """
        if isym not in _VALID_ISYM:
            raise ValueError(f"isym must be 0–4, got {isym!r}")

        header = ldt.header
        intensities = ldt.intensities

        mc = int(getattr(header, "mc", 0) or 0)
        ng = int(getattr(header, "ng", 0) or 0)
        if mc <= 0 or ng <= 0:
            raise ValueError("Ldt has empty angular geometry (mc=0 or ng=0).")

        current_isym = int(getattr(header, "isym", 0) or 0)

        # Already symmetrised and force not requested → identity copy
        if current_isym != 0 and not force:
            new_header = replace(header)
            new_intensities = [row[:] for row in intensities]
            return _rebuild_ldt(ldt, new_header, new_intensities)

        mode = _ISYM_TO_MODE[isym]

        if mode == "none":
            new_header = replace(header)
            new_intensities = [row[:] for row in intensities]
            return _rebuild_ldt(ldt, new_header, new_intensities)

        new_intensities = _apply_symmetry(header, intensities, mode, mc, ng)
        new_header = replace(header, isym=isym)
        return _rebuild_ldt(ldt, new_header, new_intensities)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _rebuild_ldt(ldt: Any, new_header: Any, new_intensities: List[List[float]]) -> Any:
    """Return a new Ldt-like object with updated header and intensities.

    Works with pyldt's Ldt dataclass (replace) or any object that exposes
    .header and .intensities attributes.
    """
    try:
        return replace(ldt, header=new_header, intensities=new_intensities)
    except TypeError:
        # Fallback for non-dataclass Ldt objects
        obj = object.__new__(type(ldt))
        obj.__dict__.update(ldt.__dict__)
        obj.header = new_header
        obj.intensities = new_intensities
        return obj


def _avg2(a: List[float], b: List[float]) -> List[float]:
    return [(x + y) * 0.5 for x, y in zip(a, b)]


def _apply_symmetry(
    header: Any,
    intensities: List[List[float]],
    mode: str,
    mc: int,
    ng: int,
) -> List[List[float]]:
    """Core averaging logic. Returns a new [mc × ng] matrix."""

    mat = [row[:] for row in intensities]
    new = [[0.0] * ng for _ in range(mc)]

    # --- ISYM 1: rotational ---
    if mode == "rot":
        base = [
            sum(mat[i][j] for i in range(mc)) / mc
            for j in range(ng)
        ]
        for i in range(mc):
            new[i] = base[:]
        return new

    # --- ISYM 2 / 3: one mirror plane ---
    if mode in ("C0", "C90"):
        Cdeg, _ = safe_angles_deg(header)
        step = 360.0 / mc if mc else 0.0
        angle_to_idx = {round(a % 360.0, 6): i for i, a in enumerate(Cdeg)}

        def idx_for(a: float) -> int:
            key = round(a % 360.0, 6)
            if key in angle_to_idx:
                return angle_to_idx[key]
            return int(round((a % 360.0) / step)) % mc

        visited = [False] * mc
        for i in range(mc):
            if visited[i]:
                continue
            C = float(Cdeg[i])
            mirror = (360.0 - C) if mode == "C0" else (180.0 - C)
            j = idx_for(mirror)
            if i == j or visited[j]:
                new[i] = mat[i][:]
                visited[i] = True
            else:
                avg = _avg2(mat[i], mat[j])
                new[i] = avg[:]
                new[j] = avg[:]
                visited[i] = visited[j] = True
        return new

    # --- ISYM 4: quadrant ---
    if mode == "both":
        half = mc // 2
        quarter = mc // 4

        # C=0° and C=180°
        if mc % 2 == 0:
            m0 = _avg2(mat[0], mat[half])
            new[0] = m0[:]
            new[half] = m0[:]
        else:
            new[0] = mat[0][:]

        # C=90° and C=270°
        if mc % 4 == 0:
            q = quarter
            qopp = (q + half) % mc
            mq = _avg2(mat[q], mat[qopp])
            new[q] = mq[:]
            new[qopp] = mq[:]

        # All other planes: average over all 4 quadrant mirrors
        for j in range(1, quarter):
            idxs = [j, (mc - j) % mc, (j + half) % mc, (half - j) % mc]
            avg = [
                sum(mat[idx][g] for idx in idxs) / 4.0
                for g in range(ng)
            ]
            for idx in idxs:
                new[idx] = avg[:]

        # Fill any remaining zeros (odd mc edge cases)
        for k in range(mc):
            if new[k] == [0.0] * ng:
                mirror = (mc - k) % mc
                new[k] = _avg2(mat[k], mat[mirror])
        return new

    raise ValueError(f"Unknown mode: {mode!r}")
