"""Single-coil air-core actuator optimization (family of `aircore_single`).

Free parameters: coil length Lc and axial position zc (coil center, relative to
magnet rest z_center = 0), magnet length Lm, coil r_in. r_out fixed at 7.75
(0.25 mm housing inside the 8.0 mm envelope radius), magnet D8 (r = 4.0) N52.

Archetypes explored:
  * base      -- the v0 candidate: long coil swallowing the magnet (coil -2..14)
  * pc (pole-centered) -- short coil centered on the magnet's TOP pole at rest;
                the other pole is pushed far away by a long magnet. Charge-model
                optimum: Lc ~ 2*a_mean ~ 12.5 mm.
  * st (straddle / "magnet between two half-strokes") -- magnet centered on the
                coil's bottom mouth at rest: each half-stroke lies on either
                side of the mouth, F ~ q*B_plateau while both poles stay in the
                antisymmetric region of the mouth transition.
  * mouth     -- magnet fully below the coil, top pole working the outside
                mouth gradient (classic half-in voice-coil position).

Protocol per the dual brief: z = [-2.5..2.5] (6 pts), exploit air-core
linearity (ni_scales -1/0/+1; force gain = (F(+1)-F(-1))/2 cancels the
mesh-noise offset exactly since F is linear in I; true cogging = 0).

Run (from repo root, inside `nix develop`, wrapped in a hard timeout):
  python sim/explore_aircore_single.py prescreen   # analytic, no FEMM
  python sim/explore_aircore_single.py r1|r2|r3    # FEMM rounds
  python sim/explore_aircore_single.py final       # best design, full protocol
"""

from __future__ import annotations

import math
import sys

from femm import Coil, Magnet, Model, MU0
from analyze import evaluate, summarize, F_TARGET

Z_SWEEP = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
SCALES_LIN = [-1.0, 0.0, 1.0]                    # air-core: F linear in I
SCALES_FULL = [-1.0, -0.5, 0.0, 0.5, 1.0]        # final-design protocol
STROKE_HALF = 2.5
CAPS_MM = 2.0
R_MAG = 4.0
R_OUT = 7.75


def single(name: str, Lc: float, zc: float, Lm: float,
           r_in: float = 4.75, r_out: float = R_OUT, ni: float = 300.0,
           **model_kw) -> Model:
    """Single-coil air-core model: coil center zc, length Lc; magnet rest at 0."""
    return Model(
        name=name,
        magnet=Magnet(radius=R_MAG, length=Lm, z_center=0.0),
        coils=[Coil("main", r_in, r_out, zc - Lc / 2, zc + Lc / 2, ni=ni)],
        **model_kw,
    )


def package_len(model: Model) -> float:
    """(stator axial extent) union (magnet extent at z = +/-2.5) + 2 mm caps."""
    hm = model.magnet.length / 2 + STROKE_HALF
    lo = min(min(c.z_bot for c in model.coils), -hm)
    hi = max(max(c.z_top for c in model.coils), +hm)
    return hi - lo + CAPS_MM


# ---------------------------------------------------------------- prescreen --

def bz_axis_annulus(z_mm: float, c: Coil, scale: float = 1.0) -> float:
    """On-axis Bz (T) of a uniform-J annular coil, standard closed form."""
    r1, r2 = c.r_in * 1e-3, c.r_out * 1e-3
    z1, z2 = c.z_bot * 1e-3, c.z_top * 1e-3
    z = z_mm * 1e-3
    j = scale * c.ni / ((r2 - r1) * (z2 - z1))          # A/m^2

    def term(zz):
        d = zz - z
        return d * math.log((r2 + math.hypot(r2, d)) / (r1 + math.hypot(r1, d)))

    return MU0 * j / 2 * (term(z2) - term(z1))


def gain_axis(model: Model, z: float) -> float:
    """Charge-model force gain (N per unit ni_scale) at magnet center z.

    F = q*(Bz(top) - Bz(bot)), q = Br*A/mu0.  On-axis approximation: ranking
    only (poles are r=4 disks; FEMM is the truth).
    """
    m = model.magnet
    q = m.br * math.pi * (m.radius * 1e-3) ** 2 / MU0
    top, bot = z + m.length / 2, z - m.length / 2
    return q * sum(bz_axis_annulus(top, c) - bz_axis_annulus(bot, c)
                   for c in model.coils)


def prescreen_row(model: Model):
    g = [gain_axis(model, z) for z in Z_SWEEP]
    p1 = sum(c.power(c.ni) for c in model.coils)        # W at ni_scale = 1
    gmin, gmax = min(g), max(g)
    km = gmin / math.sqrt(p1)
    return {
        "name": model.name, "gmin": gmin, "gmax": gmax,
        "droop_pct": 100 * (1 - gmin / gmax), "km_worst": km,
        "p_worst": (F_TARGET / km) ** 2 if km > 0 else float("inf"),
        "pkg": package_len(model), "g": g,
    }


def prescreen():
    print("analytic on-axis charge-model prescreen (ranking only)")
    print("calibration: v0 aircore_single FEMM gains were "
          "0.910/0.898/0.863/0.808/0.726/0.610 N/scale")
    base = single("base16m12", 16, 6, 12)
    r = prescreen_row(base)
    print(f"  prescreen for same geometry: "
          + "/".join(f"{x:.3f}" for x in r["g"]))
    rows = []
    # pole-centered family: coil center on top pole (zc = Lm/2 + dz)
    for Lm in (12, 16, 20):
        for Lc in (10, 12, 14, 16, 18):
            for dz in (-2, -1, 0, 1, 2):
                m = single(f"pc{Lc}m{Lm}d{dz:+d}", Lc, Lm / 2 + dz, Lm)
                if package_len(m) <= 35:
                    rows.append(prescreen_row(m))
    # straddle family: coil bottom mouth at magnet rest center (+ offset)
    for Lm in (16, 20):
        for Lc in (16, 18, 20):
            for off in (-2, 0, 2):
                m = single(f"st{Lc}m{Lm}o{off:+d}", Lc, Lc / 2 + off, Lm)
                if package_len(m) <= 35:
                    rows.append(prescreen_row(m))
    # mouth family: top pole at bottom coil mouth
    for Lm in (12, 16):
        for Lc in (10, 12, 16):
            m = single(f"mouth{Lc}m{Lm}", Lc, Lm / 2 + Lc / 2, Lm)
            if package_len(m) <= 35:
                rows.append(prescreen_row(m))
    rows.sort(key=lambda r: -r["km_worst"])
    hdr = f"{'name':16s} {'Km_w':>6s} {'P_w[W]':>7s} {'gmin':>6s} {'gmax':>6s} {'droop%':>6s} {'pkg':>5s}"
    print(hdr)
    for r in rows[:25]:
        print(f"{r['name']:16s} {r['km_worst']:6.3f} {r['p_worst']:7.2f} "
              f"{r['gmin']:6.3f} {r['gmax']:6.3f} {r['droop_pct']:6.1f} {r['pkg']:5.1f}")


# ------------------------------------------------------------- FEMM rounds --

def gains(rows):
    """Per-z force gain from the antisymmetric +/-1 pair (mesh noise cancels)."""
    byz = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
    out = {}
    for z, d in sorted(byz.items()):
        smax = max(d)
        out[z] = (d[smax] - d[-smax]) / (2 * smax), d.get(0.0)
    return out


def run_variants(models, scales, tag):
    from femm import run_sweep
    for m in models:
        rows = run_sweep(m, Z_SWEEP, scales, tag=tag, timeout=900)
        g = gains(rows)
        p1 = sum(c.power(c.ni) for c in m.coils)
        gmin = min(v[0] for v in g.values())
        gmax = max(v[0] for v in g.values())
        km = gmin / math.sqrt(p1)
        print(f"== {m.name} ==  pkg {package_len(m):.1f} mm, P(s=1) {p1:.2f} W")
        print("  z: " + "  ".join(f"{z:+.1f}:{v[0]:.3f}(n{v[1]:+.3f})"
                                  for z, v in g.items()))
        print(f"  gain min/max {gmin:.3f}/{gmax:.3f} N (droop "
              f"{100*(1-gmin/gmax):.1f}%), Km_worst {km:.3f} N/sqrtW, "
              f"P_worst(0.8N) {(F_TARGET/km)**2:.2f} W")
        ev = evaluate(m, rows)
        print(f"  analyze: P_worst {ev['power_worst_W']:.2f} W, tau "
              f"{ev['tau_s']*1e6:.0f} us, rise {ev['t_rise_s']*1e6:.0f} us, "
              f"|cog|max {ev['cogging_max_N']:.3f} N (air core: noise floor)")
        print()
        sys.stdout.flush()


def r1():
    # best prescreen members of each archetype + archetype extremes
    run_variants([
        single("pc14m20", 14, 10, 20),      # prescreen top: coil [3,17]
        single("pc16m20d1", 16, 11, 20),    # coil [3,19], flattest top group
        single("pc12m20", 12, 10, 20),      # coil [4,16], shortest good coil
        single("pc12m16d1", 12, 9, 16),     # cheaper/shorter magnet option
        single("st20m20", 20, 10, 20),      # true straddle: coil [0,20]
        single("mouth12m12", 12, 12, 12),   # pole at mouth, magnet below coil
    ], SCALES_LIN, "acs_r1")


def r2():
    # refine around round-1 winner (pole-centered, Lm=20, Km flat over
    # Lc 12..16): coil length + upward offset to symmetrize droop, Lm=18,
    # r_in sensitivity
    run_variants([
        single("pc13m20d1", 13, 11.0, 20),
        single("pc14m20d1", 14, 11.0, 20),
        single("pc15m20d1", 15, 11.0, 20),
        single("pc14m18d1", 14, 10.0, 18),         # Lm=18 midpoint
        single("pc14m20ri45", 14, 10.5, 20, r_in=4.5),
        single("pc14m20ri525", 14, 10.5, 20, r_in=5.25),
    ], SCALES_LIN, "acs_r2")


def r3():
    # symmetrized coil offset (zc = 10.5, between r1's dz=0 and r2's dz=+1)
    # + r_in curve midpoint + short-coil variant at the same offset
    run_variants([
        single("pc14m20c", 14, 10.5, 20),
        single("pc14m20c_ri46", 14, 10.5, 20, r_in=4.6),
        single("pc12m20c", 12, 10.5, 20),
    ], SCALES_LIN, "acs_r3")


def final():
    """Winner at full protocol + mesh/boundary checks (validation discipline).

    NB: Model.fine_box is declared but unused by Model.lua() (no-op), so the
    noise-robust numbers come from the antisymmetric +/-NI gain extraction
    (exact for air-core F linear in I); mesh_air only refines block labels.
    """
    import dataclasses
    m = best_model()
    run_variants([m], SCALES_FULL, "acs_final")
    mf = dataclasses.replace(m, name=m.name + "_fine", mesh_air=0.4)
    run_variants([mf], SCALES_FULL, "acs_final")
    mb = dataclasses.replace(m, name=m.name + "_r120", abc_radius=120.0)
    run_variants([mb], SCALES_LIN, "acs_final")


def best_model() -> Model:
    """The selected design (updated after each round).

    Pole-centered single coil: Lc=14 mm coil at zc=10.5 (0.5 mm above the
    magnet top pole at rest), D8x20 N52 magnet, r_in=4.75 (shared 0.75 mm
    guide-wall+clearance radial build). r_in=4.5 variant gains another +8%
    Km if a 0.5 mm radial guide build is accepted (see r2 pc14m20ri45).
    """
    return single("acs_best", 14, 10.5, 20, ni=300.0)


if __name__ == "__main__":
    for cmd in sys.argv[1:] or ["prescreen"]:
        globals()[cmd]()
