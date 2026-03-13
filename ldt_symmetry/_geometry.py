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
    """Project the intensity distribution onto a horizontal plane.

    Parameters
    ----------
    header:
        pyldt ``LdtHeader`` object (or any object with ``mc``, ``ng``,
        ``c_angles``, ``g_angles`` attributes).
    intensities:
        Full ``[mc × ng]`` intensity matrix in cd/klm.
    hemisphere:
        ``"floor"``  — γ in [0°, 90°)
        ``"ceiling"`` — γ in (90°, 180°]
    trim_pct:
        Fraction of the p95 weight used as a low-intensity threshold.
        Trims noise without distorting the centroid.
    theta_max_deg:
        Clip θ at this value before computing the projection radius.
        Avoids singularities near the horizontal plane.

    Returns
    -------
    pts : list of (weight, x, y)
        Weighted points in the projection plane.
    E : float
        Total weight (sum of all ``w``).
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
                theta_eq = theta
            else:
                if theta <= math.pi / 2:
                    continue
                theta_eq = math.pi - theta

            theta_eq = min(theta_eq, math.radians(theta_max_deg))
            w = abs(float(intensities[i][j])) * (abs(math.cos(theta_eq)) ** 3)
            r = abs(math.tan(theta_eq))
            raw.append((w, r * cos_phi, r * sin_phi))

    if not raw:
        return [], 0.0

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
        RMS radius (sqrt of mean eigenvalue sum).
    ``rho_c``
        Relative centroid offset: ``|mu| / R``.  Values > ~0.15 indicate
        a significantly off-centre distribution.
    ``kappa``
        Ellipse anisotropy: ``(λ1 − λ2) / (λ1 + λ2)``.
        0 = circular, 1 = fully elongated.
    ``alpha_deg``
        Principal axis angle in [0°, 90°] (angle vs. C=0° axis, folded).
    """
    if not pts:
        return dict(E=0.0, mu=(0.0, 0.0), R=0.0, rho_c=0.0, kappa=0.0, alpha_deg=0.0)

    W = sum(w for (w, _, _) in pts)
    mx = sum(w * x for (w, x, _) in pts) / W
    my = sum(w * y for (w, _, y) in pts) / W

    sxx = sum(w * (x - mx) ** 2 for (w, x, _) in pts) / W
    syy = sum(w * (y - my) ** 2 for (w, _, y) in pts) / W
    sxy = sum(w * (x - mx) * (y - my) for (w, x, y) in pts) / W

    disc = math.sqrt((sxx - syy) ** 2 + 4 * sxy ** 2)
    l1 = 0.5 * (sxx + syy + disc)
    l2 = 0.5 * (sxx + syy - disc)

    alpha = 0.5 * math.atan2(2 * sxy, sxx - syy)
    R = math.sqrt(l1 + l2)
    rho_c = math.sqrt(mx ** 2 + my ** 2) / (R + 1e-12)
    kappa = (l1 - l2) / (l1 + l2 + 1e-12)

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

    H4 detects 4-fold (quadrant) symmetry.  H2 detects 2-fold symmetry.
    Both are normalised to [0, 1]; higher values mean stronger anisotropy
    at the corresponding order.

    Parameters
    ----------
    pts:
        Weighted points ``(w, x, y)``.
    p:
        Radial weight exponent.  Higher values focus the analysis on the
        outer ring of the distribution.
    rband:
        ``(lo, hi)`` quantiles defining the radial band to analyse.
        (0.75, 0.98) focuses on the outer 23 % of the distribution,
        where orientation anisotropy is most visible.
    """
    if not pts:
        return {"H4": 0.0, "H2": 0.0, "phi4": 0.0, "phi2": 0.0}

    # Build (weight, radius) list for quantile computation
    wr = [(w, math.sqrt(x * x + y * y)) for (w, x, y) in pts if x * x + y * y > 0]
    if not wr:
        return {"H4": 0.0, "H2": 0.0, "phi4": 0.0, "phi2": 0.0}

    def wpercentile(data: List[Tuple[float, float]], q: float) -> float:
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

    a2c, a2s = A2c / W, A2s / W
    a4c, a4s = A4c / W, A4s / W

    return {
        "H4": math.sqrt(a4c ** 2 + a4s ** 2),
        "H2": math.sqrt(a2c ** 2 + a2s ** 2),
        "phi4": math.degrees(math.atan2(a4s, a4c)) / 4.0,
        "phi2": math.degrees(math.atan2(a2s, a2c)) / 2.0,
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

    # Off-centre distribution → asymmetric along one axis
    if rho_c > rho_th:
        if abs(alpha) <= th or abs(alpha - 90.0) <= th:
            return "C0" if abs(mx) >= abs(my) else "C90"
        return "none"

    # Centred distribution
    if kappa <= kappa_th:
        # Nearly circular → rotational or quadrant
        return "both" if H4 >= H4_th else "rot"

    # Centred but anisotropic (elongated)
    if abs(alpha) <= th or abs(alpha - 90.0) <= th:
        return "both"
    return "none"


# ---------------------------------------------------------------------------
# Combined floor + ceiling decision
# ---------------------------------------------------------------------------

# Combination table: table[dec_floor][dec_ceiling] → final mode
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

    Returns
    -------
    mode : str
        One of ``"rot"``, ``"both"``, ``"C0"``, ``"C90"``, ``"none"``.
    detail : dict
        Internal metrics (used by :class:`LdtAutoDetector` for the diag).
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
    REL_GATE = 1e-4
    if Ef > 0 and Ec > 0:
        if (Ec / Ef) < REL_GATE:
            dec_c = "none"
        elif (Ef / Ec) < REL_GATE:
            dec_f = "none"

    if Ef > 0 and Ec == 0:
        mode = dec_f
    elif Ec > 0 and Ef == 0:
        mode = dec_c
    else:
        mode = _COMBINE_TABLE.get(dec_f, {}).get(dec_c, "none")

    # Priority rule: non-circular luminous area → prefer "both" over "rot"
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

    Lower score means the distribution is more consistent with the
    requested symmetry.  A score of 0.0 means perfect symmetry.

    Parameters
    ----------
    isym:
        1 = rotational, 2 = C0–C180, 3 = C90–C270, 4 = quadrant.

    Returns
    -------
    float
        Relative RMS: ``RMS(diff) / RMS(signal)``.
        Returns ``inf`` on invalid geometry.
    """
    mc = int(getattr(header, "mc", 0) or 0)
    ng = int(getattr(header, "ng", 0) or 0)
    if mc <= 0 or ng <= 0:
        return float("inf")

    Cdeg, _ = safe_angles_deg(header)
    step = 360.0 / mc
    angle_to_idx = {round(a % 360.0, 6): i for i, a in enumerate(Cdeg)}

    def idx_for(a: float) -> int:
        key = round(a % 360.0, 6)
        if key in angle_to_idx:
            return angle_to_idx[key]
        return int(round((a % 360.0) / step)) % mc

    # ISYM 1: rotational — compare each C-plane against the mean
    if isym == 1:
        diffs: List[float] = []
        sig: List[float] = []
        for g in range(ng):
            col = [float(intensities[c][g]) for c in range(mc)]
            mu = sum(col) / mc
            diffs.extend(v - mu for v in col)
            sig.append(mu)
        return _rms(diffs) / (_rms(sig) or 1.0)

    # ISYM 2: C0–C180 — compare C with 360°−C
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

    # ISYM 3: C90–C270 — compare C with 180°−C
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

    # ISYM 4: quadrant — worst of C0–C180 and C90–C270
    if isym == 4:
        return max(
            compute_symmetry_score(header, intensities, 2),
            compute_symmetry_score(header, intensities, 3),
        )

    raise ValueError(f"isym must be 1, 2, 3 or 4 — got {isym!r}")
