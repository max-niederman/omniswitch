"""Spec-compliance evaluation of actuator sweep results.

Input: rows from femm.run_sweep (z, ni_scale, Fz, flux_<coil>...).
Output: metrics vs the hard requirements (see CLAUDE.md):
  - worst-case ohmic power to hold 0.8 N anywhere in the stroke (<= 5 W)
  - cogging force profile at zero current (steel variants)
  - bidirectionality: excitation needed for +0.8 N and -0.8 N
  - electrical time constant tau = L/R (turn-count invariant) and rise time
    to the 0.8 N current from a 24 V bus at a chosen PWM headroom
"""

from __future__ import annotations

import math

import numpy as np

from femm import Coil, Model, RHO_CU_HOT

V_BUS = 24.0
F_TARGET = 0.8
P_MAX = 5.0
SLEW_MAX_S = 10e-3
HOLD_DUTY = 0.3      # steady-state voltage fraction of bus at 0.8 N (headroom for slew)


def _by_z(rows):
    out = {}
    for r in rows:
        out.setdefault(r["z"], []).append(r)
    for z in out:
        out[z].sort(key=lambda r: r["ni_scale"])
    return out


def scale_for_force(rows_at_z, f_target):
    """ni_scale s such that Fz(s) = f_target.

    Interpolates within the swept range; linearly extrapolates from the edge
    pair outside it (exact for air-core, approximate for saturating steel —
    widen the ni_scale sweep if a steel variant lands out of range).
    """
    s = np.array([r["ni_scale"] for r in rows_at_z])
    f = np.array([r["Fz"] for r in rows_at_z])
    order = np.argsort(f)
    s, f = s[order], f[order]
    if f[0] <= f_target <= f[-1]:
        return float(np.interp(f_target, f, s))
    i, j = (0, 1) if f_target < f[0] else (-2, -1)
    slope = (s[j] - s[i]) / (f[j] - f[i])
    return float(s[j] + slope * (f_target - f[j]))


def total_power(model: Model, s: float) -> float:
    return sum(c.power(c.ni * s) for c in model.coils)


def coil_tau(model: Model, rows) -> float:
    """Electrical time constant L/R (s), turn-count invariant.

    L1 (single-turn) from d(flux)/d(NI) using the two largest ni_scale points
    at the stroke-center position; R1 (single-turn) analytic. Coils in series:
    tau = sum(L1_c) / sum(R1_c) including mutual terms is approximated by
    per-coil self terms only (conservative to ~2x; refine if marginal).
    """
    zs = sorted({r["z"] for r in rows})
    zmid = zs[len(zs) // 2]
    at_z = [r for r in rows if r["z"] == zmid]
    at_z.sort(key=lambda r: r["ni_scale"])
    hi, lo = at_z[-1], at_z[-2]
    l1_tot, r1_tot = 0.0, 0.0
    for c in model.coils:
        dni = (hi["ni_scale"] - lo["ni_scale"]) * c.ni
        if abs(dni) < 1e-12:
            continue
        dlam = hi[f"flux_{c.name}"] - lo[f"flux_{c.name}"]
        l1_tot += abs(dlam / dni)
        r1_tot += RHO_CU_HOT * 2 * math.pi * (c.r_mean_mm * 1e-3) / (
            c.area_mm2 * 1e-6 * c.fill)
    return l1_tot / r1_tot if r1_tot else 0.0


def evaluate(model: Model, rows, f_target: float = F_TARGET) -> dict:
    byz = _by_z(rows)
    zs = sorted(byz)

    cog = {z: next((r["Fz"] for r in byz[z] if r["ni_scale"] == 0), None)
           for z in zs}
    s_up = {z: scale_for_force(byz[z], +f_target) for z in zs}
    s_dn = {z: scale_for_force(byz[z], -f_target) for z in zs}
    p_up = {z: total_power(model, s_up[z]) for z in zs}
    p_dn = {z: total_power(model, s_dn[z]) for z in zs}

    p_worst = max(max(p_up.values()), max(p_dn.values()))
    tau = coil_tau(model, rows)
    # current rise to hold-point at HOLD_DUTY of bus: i/i_inf = HOLD_DUTY
    t_rise = tau * math.log(1 / (1 - HOLD_DUTY))

    return {
        "model": model.name,
        "z_positions": zs,
        "cogging_N": cog,
        "scale_up": s_up,
        "scale_down": s_dn,
        "power_up_W": p_up,
        "power_down_W": p_dn,
        "power_worst_W": p_worst,
        "km_worst": f_target / math.sqrt(p_worst) if p_worst > 0 else float("inf"),
        "tau_s": tau,
        "t_rise_s": t_rise,
        "pass_power": p_worst <= P_MAX,
        "pass_slew": t_rise <= SLEW_MAX_S,  # air-core; steel needs AC eddy check too
        "cogging_max_N": max(abs(v) for v in cog.values() if v is not None)
                         if any(v is not None for v in cog.values()) else 0.0,
    }


def summarize(ev: dict) -> str:
    zs = ev["z_positions"]
    lines = [
        f"== {ev['model']} ==",
        f"worst-case power for ±{F_TARGET} N over stroke: {ev['power_worst_W']:.2f} W "
        f"({'PASS' if ev['pass_power'] else 'FAIL'} <= {P_MAX} W)",
        f"motor constant (worst): {ev['km_worst']:.3f} N/sqrt(W)",
        f"tau = {ev['tau_s']*1e6:.0f} us, rise to hold current at {HOLD_DUTY:.0%} duty: "
        f"{ev['t_rise_s']*1e6:.0f} us ({'PASS' if ev['pass_slew'] else 'FAIL'} < 10 ms)",
        f"max |cogging|: {ev['cogging_max_N']:.3f} N",
        "  z[mm]   F=+0.8N: s, P[W]   F=-0.8N: s, P[W]   cog[N]",
    ]
    for z in zs:
        c = ev["cogging_N"][z]
        lines.append(
            f"  {z:6.2f}  {ev['scale_up'][z]:7.3f} {ev['power_up_W'][z]:6.2f}   "
            f"{ev['scale_down'][z]:8.3f} {ev['power_down_W'][z]:6.2f}   "
            f"{c if c is not None else float('nan'):7.3f}")
    return "\n".join(lines)
