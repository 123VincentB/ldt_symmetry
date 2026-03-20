# -*- coding: utf-8 -*-
"""ldt_symmetry._geometry
========================

Internal module — not part of the public API.

Provides the low-level functions used by LdtAutoDetector:
  - Angular grid reconstruction
  - Intensity projection onto floor / ceiling planes
  - Weighted moment metrics (center, anisotropy, orientation)
  - Angular harmonic decomposition (H2, H4)
  - Per-mode symmetry scores (relative RMS)
  - Shape-based mode decision (one hemisphere, then combined)

Physical background
-------------------
The auto-detection algorithm projects the luminous intensity distribution
I(C, γ) onto a horizontal plane (floor for γ < 90°, ceiling for γ > 90°)
using the standard illuminance formula for a point source::

    E = I · cos³(γ) / d²

where d is the distance to the lit surface.  The resulting illuminance
pattern (the "light patch") is then analysed geometrically:

- A circular patch suggests rotational symmetry (ISYM 1).
- A square or rectangular patch suggests quadrant symmetry (ISYM 4).
- An elongated patch aligned with C0 or C90 suggests a single mirror
  plane (ISYM 2 or 3).
- An off-centre patch reveals asymmetry along a specific axis.

Floor (direct) and ceiling (indirect) components are evaluated
independently, then combined.  Many luminaires emit only in one
hemisphere; when only one hemisphere is present the other is ignored.

The shape-based candidate is always validated by a score veto
(relative RMS of mirror-plane differences) before acceptance.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple


# ---------------------------------------------------------------------------
# Angular grid helpers
# ---------------------------------------------------------------------------

def safe_angles_deg(
    header: Any,
) -> Tuple[List[float], List[float]]:
    """Return coherent C-angle and G-angle lists.

    Uses ``header.c_angles`` / ``header.g_angles`` when available and
    consistent with ``header.mc`` / ``header.ng``.  Falls back to uniform
    grids otherwise.
    """
    mc = int(getattr(header, "mc", 0) or 0)
    ng = int(getattr(header, "ng", 0) or 0)

    c_raw = getattr(header, "c_angles", None) or []
    if c_raw and len(c_raw) == mc:
        C = list(c_raw)
    else:
        step = 360.0 / mc if mc else 0.0
        C = [i * step for i in range(mc)]

    g_raw = getattr(header, "g_angles", None) or []
    if g_raw and len(g_raw) == ng:
        G = list(g_raw)
    else:
        step = 180.0 / (ng - 1) if ng > 1 else 0.0
        G = [i * step for i in range(ng)]

    return C, G


# ---------------------------------------------------------------------------
# Intensity projection onto a plane
# ---------------------------------------------------------------------------

def plan_points(
    header: Any,
    intensities: List[List[float]],
    hemisphere: str,
    trim_pct: float = 0.02,
    theta_max_deg: float = 85.0,
) -> Tuple[List[Tuple[float, float, float]], float]:
    """Project the luminous intensity distribution onto a horizontal plane.

    Each measurement point I(C, γ) is mapped to an illuminance weight and
    a (x, y) position in the projection plane using the point-source
    illuminance formula::

        E = I · cos³(γ) / d²   →   w = I · |cos³(γ_eq)|

    where γ_eq is the polar angle measured from the vertical axis (nadir
    for floor, zenith for ceiling).  The distance d is factored out because
    only the relative shape of the light patch matters here, not its
    absolute magnitude.

    The projected position is::

        r = tan(γ_eq)          (normalised radius on the plane)
        x = r · cos(C)
        y = r · sin(C)

    This gives the shape of the illuminance footprint on the floor (or
    ceiling) as seen from directly above.  Analysing this footprint reveals
    the symmetry of the luminaire:

    - A circular footprint → rotational symmetry (ISYM 1).
    - A square / rectangular footprint → quadrant symmetry (ISYM 4).
    - An elongated footprint aligned with C0 or C90 → single mirror plane
      (ISYM 2 or 3).
    - An off-centre footprint → asymmetry along a specific axis.

    Parameters
    ----------
    header:
        pyldt ``LdtHeader`` object (or any object with ``mc``, ``ng``,
        ``c_angles``, ``g_angles`` attributes).
    intensities:
        Full ``[mc × ng]`` intensity matrix in cd/klm.
    hemisphere:
        ``"floor"``   — γ in [0°, 90°)  — direct light component.
        ``"ceiling"`` — γ in (90°, 180°] — indirect light component.
    trim_pct:
        Fraction of the p95 weight used as a low-intensity threshold.
        Trims noise at the edge of the measurement without distorting
        the centroid of the main light patch.
    theta_max_deg:
        Clip γ_eq at this value before computing tan(γ_eq).  Avoids
        singularities near the horizontal plane (γ → 90°), where the
        projection radius diverges.

    Returns
    -------
    pts : list of (weight, x, y)
        Weighted points in the projection plane.
    E : float
        Total weight (sum of all ``w``), proportional to total illuminance.
    """
    assert hemisphere in ("floor", "ceiling")

    Cdeg, Gdeg = safe_angles_deg(header)
    mc = int(getattr(header, "mc", 0) or 0)
    ng = int(getattr(header, "ng", 0) or 0)

    raw: List[Tuple[float, float, float]] = []

    for i in range(mc):
        phi = math.radians(Cdeg[i])
        cos_phi = math.cos(phi)
        sin_phi = math.sin(phi)

        for j in range(ng):
            theta = math.radians(Gdeg[j])

            if hemisphere == "floor":
                if theta > math.pi / 2:
                    continue
                # γ_eq = γ for the floor hemisphere (measured from nadir)
                theta_eq = theta
            else:
                if theta <= math.pi / 2:
                    continue
                # γ_eq = 180° − γ for the ceiling hemisphere
                # (fold back to [0°, 90°] measured from zenith)
                theta_eq = math.pi - theta

            # Clip to avoid tan() singularity near the horizontal plane
            theta_eq = min(theta_eq, math.radians(theta_max_deg))

            # Illuminance weight: I · cos³(γ_eq)  (d² factored out)
            w = abs(float(intensities[i][j])) * (abs(math.cos(theta_eq)) ** 3)

            # Projection radius: r = tan(γ_eq)  (normalised to unit distance)
            r = abs(math.tan(theta_eq))
            raw.append((w, r * cos_phi, r * sin_phi))

    if not raw:
        return [], 0.0

    # Trim low-intensity noise: discard points below trim_pct × p95
    weights_sorted = sorted(w for (w, _, _) in raw)
    p95 = weights_sorted[int(0.95 * (len(weights_sorted) - 1))]
    thr = max(0.0, trim_pct) * (p95 if p95 > 0 else 1.0)

    pts = [(w, x, y) for (w, x, y) in raw if w >= thr]
    E = sum(w for (w, _, _) in pts)
    return pts, E


# ---------------------------------------------------------------------------
# Weighted moment metrics
# ---------------------------------------------------------------------------

def moments_metrics(
    pts: List[Tuple[float, float, float]],
) -> Dict[str, Any]:
    """Compute centroid, ellipse anisotropy and principal axis orientation.

    Fits a weighted ellipse to the illuminance footprint returned by
    :func:`plan_points`.  The ellipse parameters are used to characterise
    the shape of the light patch:

    - A nearly circular ellipse (kappa ≈ 0) with a centred centroid
      (rho_c ≈ 0) suggests rotational symmetry.
    - An elongated ellipse (kappa > kappa_th) aligned with C0 or C90
      suggests a single mirror plane or quadrant symmetry.
    - A centred but circular ellipse with strong H4 harmonics (detected
      by :func:`angular_harmonics`) suggests quadrant symmetry — a square
      footprint looks circular to the ellipse fit.

    Parameters
    ----------
    pts:
        List of ``(weight, x, y)`` as returned by :func:`plan_points`.

    Returns
    -------
    dict with keys:

    ``E``
        Total weight.
    ``mu``
        Weighted centroid ``(mx, my)``.
    ``R``
        RMS radius (sqrt of the sum of both eigenvalues of the covariance
        matrix); a scale measure for the footprint size.
    ``rho_c``
        Relative centroid offset: ``|mu| / R``.  Large values (> rho_th)
        indicate that the light patch centre is shifted away from the
        goniophotometer axis, revealing asymmetry.
    ``kappa``
        Ellipse anisotropy: ``(λ1 − λ2) / (λ1 + λ2)``.
        0 = circular footprint (candidate: rot or both),
        1 = fully elongated (candidate: C0 or C90).
    ``alpha_deg``
        Principal axis angle in [0°, 90°] (folded; angle vs. C=0° axis).
        0° or 90° means the elongation is aligned with C0 or C90.
    """
    if not pts:
        return dict(E=0.0, mu=(0.0, 0.0), R=0.0, rho_c=0.0, kappa=0.0, alpha_deg=0.0)

    W = sum(w for (w, _, _) in pts)
    mx = sum(w * x for (w, x, _) in pts) / W
    my = sum(w * y for (w, _, y) in pts) / W

    # Weighted covariance matrix elements
    sxx = sum(w * (x - mx) ** 2 for (w, x, _) in pts) / W
    syy = sum(w * (y - my) ** 2 for (w, _, y) in pts) / W
    sxy = sum(w * (x - mx) * (y - my) for (w, x, y) in pts) / W

    # Eigenvalues of the 2×2 covariance matrix (analytic formula)
    disc = math.sqrt((sxx - syy) ** 2 + 4 * sxy ** 2)
    l1 = 0.5 * (sxx + syy + disc)   # major semi-axis² (largest spread)
    l2 = 0.5 * (sxx + syy - disc)   # minor semi-axis²

    # Principal axis angle (atan2 of the dominant eigenvector)
    alpha = 0.5 * math.atan2(2 * sxy, sxx - syy)

    R = math.sqrt(l1 + l2)
    rho_c = math.sqrt(mx ** 2 + my ** 2) / (R + 1e-12)
    kappa = (l1 - l2) / (l1 + l2 + 1e-12)

    # Fold alpha into [0°, 90°] — only the axis orientation matters, not
    # its direction, so C0 (0°) and C180 (180°) are equivalent here.
    alpha_deg = abs(math.degrees(alpha)) % 180.0
    if alpha_deg > 90.0:
        alpha_deg = 180.0 - alpha_deg

    return dict(E=W, mu=(mx, my), R=R, rho_c=rho_c, kappa=kappa, alpha_deg=alpha_deg)


# ---------------------------------------------------------------------------
# Angular harmonic decomposition
# ---------------------------------------------------------------------------

def angular_harmonics(
    pts: List[Tuple[float, float, float]],
    p: float = 2.0,
    rband: Tuple[float, float] = (0.0, 1.0),
) -> Dict[str, float]:
    """Compute H2 and H4 angular harmonics on a radial band.

    The ellipse fit in :func:`moments_metrics` cannot distinguish a
    circular footprint from a square one — both look circular in terms of
    eigenvalues.  Angular harmonics resolve this ambiguity:

    - H4 (4th-order harmonic) is sensitive to 4-fold (quadrant) symmetry.
      A square footprint produces a strong H4 signal even if its ellipse
      is nearly circular.
    - H2 (2nd-order harmonic) captures 2-fold symmetry (elongation),
      complementing the ellipse anisotropy kappa.

    Both are computed on the *outer radial band* of the footprint
    (default: 75th–98th percentile of the radius distribution).  The
    centre of the patch is excluded because it tends to be circular by
    nature regardless of the luminaire shape, and would dilute the signal.
    The band boundaries were tuned on the 42-sample calibration dataset.

    Both harmonics are normalised to [0, 1]; higher values indicate
    stronger anisotropy at the corresponding angular order.

    Parameters
    ----------
    pts:
        Weighted points ``(w, x, y)``.
    p:
        Radial weight exponent.  Higher values focus the analysis on the
        outer ring of the distribution.  Tuned on the 42-sample dataset.
    rband:
        ``(lo, hi)`` quantiles defining the radial band to analyse.
        Default ``(0.75, 0.98)`` focuses on the outer ~23 % of the
        footprint, where orientation anisotropy is most pronounced.
    """
    if not pts:
        return {"H4": 0.0, "H2": 0.0, "phi4": 0.0, "phi2": 0.0}

    # Build (weight, radius) list for quantile computation
    wr = [(w, math.sqrt(x * x + y * y)) for (w, x, y) in pts if x * x + y * y > 0]
    if not wr:
        return {"H4": 0.0, "H2": 0.0, "phi4": 0.0, "phi2": 0.0}

    def wpercentile(data: List[Tuple[float, float]], q: float) -> float:
        """Weighted percentile of the radius distribution."""
        total = sum(w for (w, _) in data)
        if total <= 0:
            return 0.0
        target = q * total
        acc = 0.0
        for w, r in sorted(data, key=lambda t: t[1]):
            acc += w
            if acc >= target:
                return r
        return data[-1][1]

    lo, hi = rband
    r_lo = wpercentile(wr, max(0.0, min(1.0, lo)))
    r_hi = wpercentile(wr, max(0.0, min(1.0, hi)))
    if r_hi <= r_lo:
        r_lo, r_hi = 0.0, float("inf")

    # Accumulate Fourier coefficients for orders 2 and 4
    # Each point contributes w · r^p · exp(i·n·φ) to the nth harmonic.
    # r^p up-weights the outer ring, making the analysis more sensitive
    # to the angular shape of the footprint boundary.
    A2c = A2s = A4c = A4s = W = 0.0
    for w, x, y in pts:
        r2 = x * x + y * y
        if r2 <= 0.0:
            continue
        r = math.sqrt(r2)
        if not (r_lo <= r <= r_hi):
            continue
        phi = math.atan2(y, x)
        wp = w * (r ** p)
        A2c += wp * math.cos(2.0 * phi)
        A2s += wp * math.sin(2.0 * phi)
        A4c += wp * math.cos(4.0 * phi)
        A4s += wp * math.sin(4.0 * phi)
        W += wp

    if W <= 0.0:
        return {"H4": 0.0, "H2": 0.0, "phi4": 0.0, "phi2": 0.0}

    # Normalised harmonic amplitudes (magnitude of complex Fourier coefficient)
    a2c, a2s = A2c / W, A2s / W
    a4c, a4s = A4c / W, A4s / W

    return {
        "H4": math.sqrt(a4c ** 2 + a4s ** 2),   # quadrant anisotropy (0–1)
        "H2": math.sqrt(a2c ** 2 + a2s ** 2),   # 2-fold anisotropy (0–1)
        "phi4": math.degrees(math.atan2(a4s, a4c)) / 4.0,  # H4 axis angle
        "phi2": math.degrees(math.atan2(a2s, a2c)) / 2.0,  # H2 axis angle
    }


# ---------------------------------------------------------------------------
# Per-hemisphere mode decision
# ---------------------------------------------------------------------------

def _decide_hemisphere(
    metrics: Dict[str, Any],
    rho_th: float,
    kappa_th: float,
    theta_ax_deg: float,
    H4_th: float,
) -> str:
    """Return the best symmetry mode for one hemisphere.

    Decision logic (applied to the illuminance footprint metrics):

    1. **Off-centre footprint** (rho_c > rho_th): the light patch centroid
       is shifted away from the goniophotometer axis.  The shift direction
       determines the axis of asymmetry → C0 or C90 (single mirror plane).
       If the shift is diagonal (not aligned with C0 or C90) → ``"none"``.

    2. **Centred, circular footprint** (rho_c ≤ rho_th, kappa ≤ kappa_th):
       the patch is nearly circular.  H4 distinguishes a true circle
       (rotational symmetry) from a square footprint (quadrant symmetry):
       - H4 ≥ H4_th → ``"both"`` (square/rectangular footprint)
       - H4 < H4_th → ``"rot"``  (circular footprint)

    3. **Centred, elongated footprint** (rho_c ≤ rho_th, kappa > kappa_th):
       the patch is elliptical.  If the major axis is aligned with C0 or
       C90 (within theta_ax_deg) → ``"both"`` (rectangular footprint
       with axis-aligned elongation).  Otherwise → ``"none"``.

    The kappa_th threshold (default 0.05) is intentionally low: even a
    slight flattening of the footprint is enough to rule out pure
    rotational symmetry.  Tuned on the 42-sample calibration dataset.

    Returns one of: ``"rot"``, ``"both"``, ``"C0"``, ``"C90"``, ``"none"``.
    """
    if metrics.get("E", 0.0) <= 0.0:
        return "none"

    rho_c = float(metrics["rho_c"])
    kappa = float(metrics["kappa"])
    alpha = float(metrics["alpha_deg"])
    mx, my = metrics["mu"]
    H4 = float(metrics.get("H4", 0.0))
    th = float(theta_ax_deg)

    # --- Case 1: off-centre footprint → asymmetry along one axis ---
    if rho_c > rho_th:
        if abs(alpha) <= th or abs(alpha - 90.0) <= th:
            # Shift is axis-aligned: the larger shift component determines
            # which mirror plane is relevant (C0 ↔ x-axis, C90 ↔ y-axis)
            return "C0" if abs(mx) >= abs(my) else "C90"
        # Diagonal shift: no simple mirror plane can be identified
        return "none"

    # --- Case 2: centred, circular footprint ---
    if kappa <= kappa_th:
        # H4 distinguishes a circular footprint (rot) from a square one (both)
        return "both" if H4 >= H4_th else "rot"

    # --- Case 3: centred, elongated footprint ---
    # An ellipse aligned with C0 (alpha ≈ 0°) or C90 (alpha ≈ 90°) is
    # consistent with a rectangular luminous area → quadrant symmetry.
    if abs(alpha) <= th or abs(alpha - 90.0) <= th:
        return "both"
    # Diagonal elongation: not attributable to a standard symmetry mode
    return "none"


# ---------------------------------------------------------------------------
# Combined floor + ceiling decision
# ---------------------------------------------------------------------------

# Combination table: table[dec_floor][dec_ceiling] → final mode
#
# A luminaire may emit light in one or both hemispheres:
#   - floor (direct):   γ ∈ [0°, 90°)  — downward light
#   - ceiling (indirect): γ ∈ (90°, 180°] — upward light
#
# Each hemisphere is evaluated independently and produces a candidate mode.
# The two candidates are then combined using the following rules:
#
# - When both hemispheres agree (diagonal) → that mode wins.
# - When one hemisphere is absent ("none") → the other hemisphere's mode
#   prevails.  Example: a purely direct luminaire has no ceiling component;
#   the floor decision is used as-is.
# - When the two hemispheres disagree → the more restrictive (less
#   symmetric) mode wins.  Rationale: if the floor sees C0 symmetry but
#   the ceiling sees full rotation, the C0 constraint from the floor must
#   be respected.  Examples:
#     rot  + C0  → C0    (C0 is more restrictive than rot)
#     both + C0  → C0    (C0 is more restrictive than both)
#     C0   + C90 → both  (two orthogonal mirror planes → quadrant)
_COMBINE_TABLE = {
    "rot":  {"rot": "rot",  "both": "both", "C0": "C0",   "C90": "C90",  "none": "rot"},
    "both": {"rot": "both", "both": "both", "C0": "C0",   "C90": "C90",  "none": "both"},
    "C0":   {"rot": "C0",   "both": "C0",   "C0": "C0",   "C90": "both", "none": "C0"},
    "C90":  {"rot": "C90",  "both": "C90",  "C0": "both", "C90": "C90",  "none": "C90"},
    "none": {"rot": "rot",  "both": "both", "C0": "C0",   "C90": "C90",  "none": "none"},
}


def decide_by_shape(
    header: Any,
    intensities: List[List[float]],
    rho_th: float = 0.04,
    kappa_th: float = 0.05,
    theta_ax_deg: float = 15.0,
    trim_pct: float = 0.02,
    theta_max_deg: float = 85.0,
    H4_th: float = 0.045,
    H4_p: float = 5.0,
    H4_rband: Tuple[float, float] = (0.75, 0.98),
) -> Tuple[str, Dict[str, Any]]:
    """Shape-based symmetry mode decision (floor + ceiling combined).

    Projects I(C, γ) onto floor and ceiling planes, analyses each
    illuminance footprint independently, then combines the two decisions
    using :data:`_COMBINE_TABLE`.

    An energy gate suppresses a hemisphere whose total illuminance weight
    is negligible relative to the other (ratio < 1e-4).  This avoids
    noise in the near-zero hemisphere biasing the decision.

    Returns
    -------
    mode : str
        One of ``"rot"``, ``"both"``, ``"C0"``, ``"C90"``, ``"none"``.
    detail : dict
        Internal metrics for both hemispheres (used by
        :class:`LdtAutoDetector` for the ``diag`` output).
    """
    pts_f, Ef = plan_points(header, intensities, "floor",
                            trim_pct=trim_pct, theta_max_deg=theta_max_deg)
    met_f = moments_metrics(pts_f)
    met_f.update(angular_harmonics(pts_f, p=H4_p, rband=H4_rband))
    dec_f = _decide_hemisphere(met_f, rho_th, kappa_th, theta_ax_deg, H4_th)

    pts_c, Ec = plan_points(header, intensities, "ceiling",
                            trim_pct=trim_pct, theta_max_deg=theta_max_deg)
    met_c = moments_metrics(pts_c)
    met_c.update(angular_harmonics(pts_c, p=H4_p, rband=H4_rband))
    dec_c = _decide_hemisphere(met_c, rho_th, kappa_th, theta_ax_deg, H4_th)

    # Energy gate: ignore a hemisphere whose energy is negligible
    # (ratio < 1e-4 of the dominant hemisphere)
    REL_GATE = 1e-4
    if Ef > 0 and Ec > 0:
        if (Ec / Ef) < REL_GATE:
            dec_c = "none"   # ceiling energy negligible → ignore
        elif (Ef / Ec) < REL_GATE:
            dec_f = "none"   # floor energy negligible → ignore

    if Ef > 0 and Ec == 0:
        mode = dec_f
    elif Ec > 0 and Ef == 0:
        mode = dec_c
    else:
        mode = _COMBINE_TABLE.get(dec_f, {}).get(dec_c, "none")

    # Priority rule: a non-circular luminous area (width_lum_area > 0)
    # cannot produce a truly rotationally symmetric footprint — the
    # physical geometry of the source imposes at least quadrant symmetry.
    # Upgrade "rot" to "both" in that case.
    width_lum = getattr(header, "width_lum_area", None)
    priority_applied = False
    if width_lum is not None and float(width_lum) > 0 and mode == "rot":
        mode = "both"
        priority_applied = True

    detail = {
        "floor": {"decision": dec_f, "E": Ef, **met_f},
        "ceiling": {"decision": dec_c, "E": Ec, **met_c},
        "priority_rule_applied": priority_applied,
    }
    return mode, detail


# ---------------------------------------------------------------------------
# Symmetry score (relative RMS)
# ---------------------------------------------------------------------------

def _rms(arr: List[float]) -> float:
    if not arr:
        return 0.0
    return math.sqrt(sum(x * x for x in arr) / len(arr))


def compute_symmetry_score(
    header: Any,
    intensities: List[List[float]],
    isym: int,
) -> float:
    """Compute a relative RMS asymmetry score for the given ISYM mode.

    Measures how well the intensity distribution actually satisfies the
    mirror (or rotational) symmetry implied by ``isym``.  Used as a veto
    in :class:`LdtAutoDetector`: if the score exceeds the threshold, the
    shape-based candidate is rejected and ISYM 0 is returned instead.

    For each symmetry mode, the score compares pairs of C-planes that
    should be identical if the symmetry were perfect:

    - **ISYM 1 (rotational)**: every C-plane is compared against the
      mean of all C-planes.  A perfectly rotationally symmetric
      distribution has zero variance across C-planes.
    - **ISYM 2 (C0–C180 mirror)**: C[i] is compared against C[360°−i].
    - **ISYM 3 (C90–C270 mirror)**: C[i] is compared against C[180°−i].
    - **ISYM 4 (quadrant)**: worst of ISYM 2 and ISYM 3 scores.  Both
      mirror planes must pass for the quadrant mode to be accepted.

    The score is the relative RMS::

        score = RMS(diff) / RMS(signal)

    where ``diff`` is the element-wise difference between mirror planes
    and ``signal`` is their mean absolute value.  A score of 0 means
    perfect symmetry; a score of 1 means the differences are as large
    as the signal itself.

    Parameters
    ----------
    isym:
        1 = rotational, 2 = C0–C180, 3 = C90–C270, 4 = quadrant.

    Returns
    -------
    float
        Relative RMS score ≥ 0.  Returns ``inf`` on invalid geometry.
    """
    mc = int(getattr(header, "mc", 0) or 0)
    ng = int(getattr(header, "ng", 0) or 0)
    if mc <= 0 or ng <= 0:
        return float("inf")

    Cdeg, _ = safe_angles_deg(header)
    step = 360.0 / mc
    angle_to_idx = {round(a % 360.0, 6): i for i, a in enumerate(Cdeg)}

    def idx_for(a: float) -> int:
        """Return the matrix column index for C-angle a (degrees)."""
        key = round(a % 360.0, 6)
        if key in angle_to_idx:
            return angle_to_idx[key]
        return int(round((a % 360.0) / step)) % mc

    # ISYM 1: rotational — compare each C-plane against the mean of all
    if isym == 1:
        diffs: List[float] = []
        sig: List[float] = []
        for g in range(ng):
            col = [float(intensities[c][g]) for c in range(mc)]
            mu = sum(col) / mc
            diffs.extend(v - mu for v in col)
            sig.append(mu)
        return _rms(diffs) / (_rms(sig) or 1.0)

    # ISYM 2: C0–C180 mirror — compare C[i] with C[360°−i]
    if isym == 2:
        visited = [False] * mc
        diffs = []
        sig = []
        for i in range(mc):
            if visited[i]:
                continue
            j = idx_for(360.0 - Cdeg[i])
            if j == i:
                visited[i] = True
                continue
            visited[i] = visited[j] = True
            for g in range(ng):
                a, b = float(intensities[i][g]), float(intensities[j][g])
                diffs.append(a - b)
                sig.append(0.5 * (abs(a) + abs(b)))
        return _rms(diffs) / (_rms(sig) or 1.0)

    # ISYM 3: C90–C270 mirror — compare C[i] with C[180°−i]
    if isym == 3:
        visited = [False] * mc
        diffs = []
        sig = []
        for i in range(mc):
            if visited[i]:
                continue
            j = idx_for(180.0 - Cdeg[i])
            if j == i:
                visited[i] = True
                continue
            visited[i] = visited[j] = True
            for g in range(ng):
                a, b = float(intensities[i][g]), float(intensities[j][g])
                diffs.append(a - b)
                sig.append(0.5 * (abs(a) + abs(b)))
        return _rms(diffs) / (_rms(sig) or 1.0)

    # ISYM 4: quadrant — both mirror planes must satisfy their threshold;
    # the worst score determines acceptance
    if isym == 4:
        return max(
            compute_symmetry_score(header, intensities, 2),
            compute_symmetry_score(header, intensities, 3),
        )

    raise ValueError(f"isym must be 1, 2, 3 or 4 — got {isym!r}")
