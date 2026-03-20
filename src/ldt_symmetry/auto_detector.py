# -*- coding: utf-8 -*-
"""ldt_symmetry.auto_detector
==============================

Provides :class:`LdtAutoDetector`, an optional helper class that
automatically detects the most appropriate ISYM symmetry mode for a
given ``Ldt`` object.

This class is intentionally separate from :class:`~ldt_symmetry.LdtSymmetriser`.
Use it only when the symmetry mode is not known in advance.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

from ._geometry import compute_symmetry_score, decide_by_shape

# Internal mode string → ISYM integer
_MODE_TO_ISYM = {"none": 0, "rot": 1, "C0": 2, "C90": 3, "both": 4}


class LdtAutoDetector:
    """Automatically detect the ISYM symmetry mode of a ``Ldt`` object.

    Detection proceeds in two stages:

    1. **Shape analysis** — projects the intensity distribution onto the
       floor and ceiling planes, computes weighted moments and angular
       harmonics, and selects a candidate mode.

    2. **Score veto** — measures the actual asymmetry for the candidate
       mode (relative RMS).  If the score exceeds the threshold the
       candidate is rejected and ISYM 0 is returned.

    The default thresholds were tuned on a dataset of real manufacturer
    LDT files and are correct in the vast majority of cases.  They can
    be adjusted via the constructor if needed.

    Parameters
    ----------
    score_th : float
        Relative RMS threshold for ISYM 2, 3, 4.
        Above this value the candidate is vetoed.  Default: 0.23.
    rot_score_th : float
        Relative RMS threshold for ISYM 1 (rotational).
        Rotational symmetry is stricter by nature.  Default: 0.20.
    rho_th : float
        Shape param — centroid offset threshold.  Default: 0.04.
    kappa_th : float
        Shape param — ellipse anisotropy threshold.  Default: 0.05.
    theta_ax_deg : float
        Shape param — axis alignment tolerance (degrees).  Default: 15.0.
    trim_pct : float
        Shape param — low-intensity trim fraction.  Default: 0.02.
    theta_max_deg : float
        Shape param — projection clip angle (degrees).  Default: 85.0.
    H4_th : float
        Shape param — H4 harmonic threshold for quadrant detection.
        Default: 0.045.
    H4_p : float
        Shape param — radial weight exponent for H4.  Default: 5.0.
    H4_rband : tuple of float
        Shape param — radial band (lo, hi quantiles) for H4.
        Default: (0.75, 0.98).

    Usage
    -----
    ::

        from pyldt import LdtReader, LdtWriter
        from ldt_symmetry import LdtAutoDetector, LdtSymmetriser

        ldt = LdtReader.read("luminaire.ldt")

        detector = LdtAutoDetector()
        isym = detector.detect(ldt)
        ldt_sym = LdtSymmetriser.symmetrise(ldt, isym=isym)
        LdtWriter.write(ldt_sym, "luminaire_SYM.ldt")

        # With diagnostic information:
        isym, diag = detector.detect(ldt, return_diag=True)
        print(diag)
    """

    def __init__(
        self,
        *,
        score_th: float = 0.23,
        rot_score_th: float = 0.20,
        rho_th: float = 0.04,
        kappa_th: float = 0.05,
        theta_ax_deg: float = 15.0,
        trim_pct: float = 0.02,
        theta_max_deg: float = 85.0,
        H4_th: float = 0.045,
        H4_p: float = 5.0,
        H4_rband: Tuple[float, float] = (0.75, 0.98),
    ) -> None:
        self.score_th = score_th
        self.rot_score_th = rot_score_th
        self.rho_th = rho_th
        self.kappa_th = kappa_th
        self.theta_ax_deg = theta_ax_deg
        self.trim_pct = trim_pct
        self.theta_max_deg = theta_max_deg
        self.H4_th = H4_th
        self.H4_p = H4_p
        self.H4_rband = H4_rband

    def detect(
        self,
        ldt: Any,
        *,
        return_diag: bool = False,
    ) -> int | Tuple[int, Dict[str, Any]]:
        """Detect the most appropriate ISYM symmetry mode.

        Parameters
        ----------
        ldt:
            A ``pyldt.Ldt`` object as returned by ``LdtReader.read()``.
        return_diag:
            If ``True``, return a ``(isym, diag)`` tuple instead of just
            ``isym``.  The ``diag`` dict contains the shape metrics, symmetry
            scores, and the veto decision.  Useful for debugging or for
            understanding why a particular mode was chosen or rejected.

        Returns
        -------
        isym : int
            Detected ISYM mode (0–4).  Returns 0 if no symmetric mode is
            plausible.
        diag : dict
            Only returned when ``return_diag=True``.  Keys:

            ``isym_detected``
                Final detected ISYM (same as the first return value).
            ``isym_candidate``
                ISYM proposed by the shape stage before veto.
            ``accepted``
                ``True`` if the candidate passed the score veto.
            ``scores``
                Dict of relative RMS scores keyed by ISYM integer
                (1, 2, 3, 4).
            ``score_threshold``
                Threshold applied to ISYM 2, 3, 4.
            ``rot_score_threshold``
                Threshold applied to ISYM 1.
            ``shape_detail``
                Internal shape metrics (floor/ceiling moments, harmonics).

        Raises
        ------
        ValueError
            If the ``Ldt`` object has empty angular geometry.
        """
        header = ldt.header
        intensities = ldt.intensities

        mc = int(getattr(header, "mc", 0) or 0)
        ng = int(getattr(header, "ng", 0) or 0)
        if mc <= 0 or ng <= 0:
            raise ValueError("Ldt has empty angular geometry (mc=0 or ng=0).")

        # Stage 1: shape-based candidate
        candidate_mode, shape_detail = decide_by_shape(
            header, intensities,
            rho_th=self.rho_th,
            kappa_th=self.kappa_th,
            theta_ax_deg=self.theta_ax_deg,
            trim_pct=self.trim_pct,
            theta_max_deg=self.theta_max_deg,
            H4_th=self.H4_th,
            H4_p=self.H4_p,
            H4_rband=self.H4_rband,
        )
        candidate_isym = _MODE_TO_ISYM[candidate_mode]

        # Stage 2: score veto
        scores = {
            1: compute_symmetry_score(header, intensities, 1),
            2: compute_symmetry_score(header, intensities, 2),
            3: compute_symmetry_score(header, intensities, 3),
            4: compute_symmetry_score(header, intensities, 4),
        }

        accepted = False
        if candidate_isym != 0:
            accepted = self._passes_veto(candidate_isym, scores)

        final_isym = candidate_isym if accepted else 0

        if not return_diag:
            return final_isym

        diag: Dict[str, Any] = {
            "isym_detected": final_isym,
            "isym_candidate": candidate_isym,
            "accepted": accepted,
            "scores": scores,
            "score_threshold": self.score_th,
            "rot_score_threshold": self.rot_score_th,
            "shape_detail": shape_detail,
        }
        return final_isym, diag

    def _passes_veto(self, isym: int, scores: Dict[int, float]) -> bool:
        """Return True if the score for ``isym`` is within the threshold."""
        if isym == 1:
            return scores[1] <= self.rot_score_th
        if isym == 2:
            return scores[2] <= self.score_th
        if isym == 3:
            return scores[3] <= self.score_th
        if isym == 4:
            # Quadrant requires both C0-C180 and C90-C270 to pass
            return scores[2] <= self.score_th and scores[3] <= self.score_th
        return False
