"""The three finalist actuator designs, as clean importable builders.

Consolidated from the exploration phase (2026-07-28); see results/EXPLORATION.md
for the full digest and sim/explore_*.py for the search history. All three
pass every hard requirement (CLAUDE.md) with N52 at 20 C:

  opp24c3      opposed pair of N52 D8x12 (like N poles facing, 0.2 mm spacer)
               mover; 3 air-core coils. 0.36 W worst hold, 35.0 mm, 0 cogging.
  w_o3_m12_r   opposed pair N52 D8x12 + 2 mm steel washer between the facing
               N poles (mover); ironless 3-coil stator. 0.42-0.48 W, 33.0 mm.
  st_w03       single N52 D8x12 mover, dual coil, 0.3 mm 1018 shell doubling
               as housing. 1.01 W worst hold, 28.0 mm, cogging ~0.07 N.

Each builder returns (model, meta): a ready-to-sweep model object and a dict
of the verified exploration numbers (worst-case hold power, package length,
base NI per coil, force-gain-vs-z table from the explore CSVs).

The model classes are inlined below (NOT imported from explore_* files, which
are frozen history) with attribution comments; parameter values are copied
EXACTLY from the explore-phase finalists.

Verification caveats that carry over (results/EXPLORATION.md "Caveats"):
  * weighted-stress-tensor noise floor ~0.03-0.06 N with automesh; use fine
    mesh and antisymmetric +/-s differencing for small forces (cogging);
  * Coil fill=0.6 is 10-20 % optimistic; quote powers at fill=0.5 too
    (meta["worst_power_fill05_W"]);
  * abc_radius=60 is converged (<=0.3 % force shift vs 120).

Run (inside `nix develop`, from repo root, ALWAYS under `timeout`):
    python sim/designs.py <name> [z1,z2,...] [s1,s2,...] [timeout_s]
    e.g. timeout 900 nix develop . -c python sim/designs.py st_w03 -2.5,0,2.5 0,1 800
Prints analyze.summarize plus a per-z force-gain comparison against the
explore-phase reference table in meta.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import MU0, Coil, Magnet, Model, Steel, run_sweep, z_path  # noqa: E402
from analyze import evaluate, summarize  # noqa: E402


# =====================================================================
# Multi-part-mover model support
# (copied from sim/explore_wildcard.py: MagnetPart / MoverSteel / MoverModel,
#  DC parts only — the AC-harmonic machinery there is dead with this
#  toolchain, FEMM harmonic solves hang headless under Wine.)
# =====================================================================

@dataclass
class MagnetPart:
    """Axially magnetized cylinder; z_offset is the part center relative to
    the mover reference z (the swept coordinate). magdir 90 = +z, 270 = -z."""
    radius: float
    length: float
    z_offset: float
    magdir: float = 90.0
    br: float = 1.43
    mu_r: float = 1.05

    @property
    def hc(self) -> float:
        return self.br / (self.mu_r * MU0)


@dataclass
class MoverSteel:
    """Steel part on the mover; z coords relative to mover reference z."""
    name: str
    r_in: float
    r_out: float
    z_bot: float
    z_top: float
    material: str = "1018 Steel"


@dataclass
class MoverModel(Model):
    """Model with a multi-part mover: list of MagnetPart + mover steels that
    translate with z (base-class mover Steel stays at fixed absolute z, which
    is wrong for washers). Also implements fine_box (declared but unused in
    the base class): an air rect labelled at mesh_air wrapped around the
    mover so the stress-tensor integration air is finely meshed. The fine box
    omits its axis edge (open_axis) to avoid overlapping collinear segments
    with magnet axis edges — identical shared edges between butted parts are
    deduped python-side."""
    magnets: list[MagnetPart] = field(default_factory=list)
    mover_steels: list[MoverSteel] = field(default_factory=list)

    def _emit_doc(self, w, z: float):
        segs: set = set()

        def rect(r1, z1, r2, z2, open_axis=False):
            edges = [(r1, z1, r2, z1), (r2, z1, r2, z2), (r2, z2, r1, z2)]
            if not open_axis:
                edges.append((r1, z2, r1, z1))
            for (a, b, c, d) in edges:
                key = tuple(sorted([(round(a, 6), round(b, 6)),
                                    (round(c, 6), round(d, 6))]))
                if key in segs:
                    continue          # butted parts share identical edges
                segs.add(key)
                w(f"mi_addnode({a:.6g}, {b:.6g})")
                w(f"mi_addnode({c:.6g}, {d:.6g})")
                w(f"mi_addsegment({a:.6g}, {b:.6g}, {c:.6g}, {d:.6g})")

        def label(r, zz, mat, circuit="", magdir=0, group=0, turns=0, mesh=0.0):
            w(f"mi_addblocklabel({r:.6g}, {zz:.6g})")
            w(f"mi_selectlabel({r:.6g}, {zz:.6g})")
            automesh = 1 if mesh <= 0 else 0
            w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, "{circuit}", '
              f"{magdir}, {group}, {turns})")
            w("mi_clearselected()")

        w("newdocument(0)")
        w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
        w('mi_getmaterial("Air")')
        w('mi_getmaterial("Copper")')
        for i, mp in enumerate(self.magnets):
            w(f'mi_addmaterial("Magnet{i}", {mp.mu_r}, {mp.mu_r}, {mp.hc:.1f}, '
              "0, 0.667, 0, 0, 1, 0, 0, 0)")
        seen = set()
        for st in list(self.mover_steels) + list(self.steels):
            if st.material not in seen:
                w(f'mi_getmaterial("{st.material}")')
                seen.add(st.material)

        for i, mp in enumerate(self.magnets):
            zc = z + mp.z_offset
            rect(0.0, zc - mp.length / 2, mp.radius, zc + mp.length / 2)
            label(mp.radius / 2, zc, f"Magnet{i}", magdir=mp.magdir, group=1,
                  mesh=self.mesh_air)

        for st in self.mover_steels:
            rect(st.r_in, z + st.z_bot, st.r_out, z + st.z_top)
            label((st.r_in + st.r_out) / 2, z + (st.z_bot + st.z_top) / 2,
                  st.material, group=1, mesh=self.mesh_air)

        for c in self.coils:
            rect(c.r_in, c.z_bot, c.r_out, c.z_top)
            w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
            label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
                  circuit=c.name, turns=1, mesh=self.mesh_air)

        for st in self.steels:   # stator-fixed steel
            rect(st.r_in, st.z_bot, st.r_out, st.z_top)
            label((st.r_in + st.r_out) / 2, (st.z_bot + st.z_top) / 2,
                  st.material, group=1 if st.mover else 0, mesh=self.mesh_air)

        if self.fine_box:
            fr, fz1, fz2 = self.fine_box
            rect(0.0, fz1, fr, fz2, open_axis=True)
            label(fr / 2, fz2 - 0.5, "Air", mesh=self.mesh_air)

        w(f"mi_makeABC(7, {self.abc_radius:.6g}, 0, 0, 0)")
        label(self.abc_radius * 0.5, self.abc_radius * 0.5, "Air")

    def lua(self, z_positions, ni_scales, out_csv, fem_file) -> str:
        L: list[str] = []
        w = L.append
        w(f'LOG = "{z_path(out_csv)}.log"')
        w("function mark(s)")
        w('    local h = openfile(LOG, "a")')
        w('    write(h, s, "\\n")')
        w("    closefile(h)")
        w("end")
        w(f'handle = openfile("{z_path(out_csv)}", "w")')
        flux_cols = "".join(f",flux_{c.name}" for c in self.coils)
        w(f'write(handle, "z,ni_scale,Fz{flux_cols}\\n")')
        for z in z_positions:
            self._emit_doc(w, z)
            w(f'mark("z={z:.6g} built")')
            w(f'mi_saveas("{z_path(fem_file)}")')
            for s in ni_scales:
                for c in self.coils:
                    w(f'mi_modifycircprop("{c.name}", 1, {c.ni * s:.6g})')
                w(f'mark("z={z:.6g} s={s:.6g} analyze")')
                w("mi_analyze()")
                w('mark("solved")')
                w("mi_loadsolution()")
                w("mo_groupselectblock(1)")
                w("fz = mo_blockintegral(19)")
                w(f'write(handle, {z:.6g}, ",", {s:.6g}, ",", fz)')
                for c in self.coils:
                    w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
                    w('write(handle, ",", lam)')
                w('write(handle, "\\n")')
                w("mo_close()")
            w("mi_close()")
        w("closefile(handle)")
        w("quit()")
        return "\n".join(L) + "\n"


# =====================================================================
# Steel-shell model with fine-mesh support
# (copied from sim/explore_steel_shell.py: ShellModel.)
# =====================================================================

@dataclass
class ShellModel(Model):
    """Model + connected steel polygons + optional fine air annulus.

    steel_polys: list of (name, [(r, z), ...], (label_r, label_z), material)
        drawn as a closed polygon, stator (group 0).
    fine_annulus: (r_in, r_out, z_min, z_max, mesh_mm) air annulus in the
        bore gap given an explicit mesh size (base Model.mesh_air only sizes
        magnet/coil/steel labels; surrounding air stays automeshed).
    """
    steel_polys: list = field(default_factory=list)
    fine_annulus: tuple | None = None

    def lua(self, z_positions, ni_scales, out_csv, fem_file) -> str:
        m = self.magnet
        L: list[str] = []
        w = L.append

        w(f'LOG = "{z_path(out_csv)}.log"')
        w("function mark(s)")
        w("    local h = openfile(LOG, \"a\")")
        w('    write(h, s, "\\n")')
        w("    closefile(h)")
        w("end")
        w(f'handle = openfile("{z_path(out_csv)}", "w")')
        flux_cols = "".join(f",flux_{c.name}" for c in self.coils)
        w(f'write(handle, "z,ni_scale,Fz{flux_cols}\\n")')

        def seg(r1, z1, r2, z2):
            w(f"mi_addnode({r1:.6g}, {z1:.6g})")
            w(f"mi_addnode({r2:.6g}, {z2:.6g})")
            w(f"mi_addsegment({r1:.6g}, {z1:.6g}, {r2:.6g}, {z2:.6g})")

        def rect(r1, z1, r2, z2):
            for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                                 (r2, z2, r1, z2), (r1, z2, r1, z1)]:
                seg(a, b, c, d)

        def poly(pts):
            n = len(pts)
            for i in range(n):
                r1, z1 = pts[i]
                r2, z2 = pts[(i + 1) % n]
                seg(r1, z1, r2, z2)

        def label(r, z, mat, circuit="", magdir=0, group=0, turns=0, mesh=0.0):
            w(f"mi_addblocklabel({r:.6g}, {z:.6g})")
            w(f"mi_selectlabel({r:.6g}, {z:.6g})")
            automesh = 1 if mesh <= 0 else 0
            w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, "{circuit}", '
              f"{magdir}, {group}, {turns})")
            w("mi_clearselected()")

        for z in z_positions:
            w("newdocument(0)")
            w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
            w('mi_getmaterial("Air")')
            w('mi_getmaterial("Copper")')
            w(f'mi_addmaterial("Magnet", {m.mu_r}, {m.mu_r}, {m.hc:.1f}, '
              "0, 0.667, 0, 0, 1, 0, 0, 0)")
            mats = {s.material for s in self.steels}
            mats |= {p[3] for p in self.steel_polys}
            for mat in sorted(mats):
                w(f'mi_getmaterial("{mat}")')

            mz1, mz2 = z - m.length / 2, z + m.length / 2
            rect(0, mz1, m.radius, mz2)
            label(m.radius / 2, z, "Magnet", magdir=90, group=1,
                  mesh=self.mesh_air)

            for c in self.coils:
                rect(c.r_in, c.z_bot, c.r_out, c.z_top)
                w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
                label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
                      circuit=c.name, turns=1, mesh=self.mesh_air)

            for s in self.steels:
                rect(s.r_in, s.z_bot, s.r_out, s.z_top)
                label((s.r_in + s.r_out) / 2, (s.z_bot + s.z_top) / 2,
                      s.material, group=1 if s.mover else 0,
                      mesh=self.mesh_air)

            for (name, pts, (lr, lz), mat) in self.steel_polys:
                poly(pts)
                label(lr, lz, mat, group=0, mesh=self.mesh_air)

            if self.fine_annulus is not None:
                fr1, fr2, fz1, fz2, fmesh = self.fine_annulus
                rect(fr1, fz1, fr2, fz2)
                label((fr1 + fr2) / 2, (fz1 + fz2) / 2, "Air", mesh=fmesh)

            w(f"mi_makeABC(7, {self.abc_radius:.6g}, 0, 0, 0)")
            label(self.abc_radius * 0.5, self.abc_radius * 0.5, "Air")

            w(f'mark("z={z:.6g} built")')
            w(f'mi_saveas("{z_path(fem_file)}")')
            for s in ni_scales:
                for c in self.coils:
                    w(f'mi_modifycircprop("{c.name}", 1, {c.ni * s:.6g})')
                w(f'mark("z={z:.6g} s={s:.6g} analyze")')
                w("mi_analyze()")
                w('mark("solved")')
                w("mi_loadsolution()")
                w("mo_groupselectblock(1)")
                w("fz = mo_blockintegral(19)")
                w(f'write(handle, {z:.6g}, ",", {s:.6g}, ",", fz)')
                for c in self.coils:
                    w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
                    w('write(handle, ",", lam)')
                w('write(handle, "\\n")')
                w("mo_close()")
            w("mi_close()")
        w("closefile(handle)")
        w("quit()")
        return "\n".join(L) + "\n"


# =====================================================================
# Finalist builders. Parameter values are copied EXACTLY from the
# explore-phase winners; gain tables are (F(+1) - F(-1))/2 per z from the
# named results/ CSVs (antisymmetric differencing cancels the constant
# stress-tensor mesh offset).
# =====================================================================

def opp24c3():
    """Finalist 1 — all-air opposed-pair, 3-coil (sim/explore_aircore_dual.py
    variant "opp24c3", stack [(-6.1, 12, +1), (+6.1, 12, -1)]).

    Two N52 D8x12, like (N) poles facing across a 0.2 mm spacer, on the
    mover; lo/mid/hi air-core coils at r 4.5-7.75 (0.5 mm guide wall +
    0.25 mm housing = OD 16.0). Originally expressed via that file's private
    Variant/gen_lua stack generator; re-expressed here as a MoverModel with
    two MagnetParts — geometrically identical (same rectangles, magdirs,
    materials).  Zero cogging (no steel anywhere), F exactly linear in NI.
    """
    model = MoverModel(
        name="opp24c3",
        magnet=Magnet(radius=4.0, length=24.2, z_center=0.0),  # placeholder, unused
        coils=[Coil("lo", 4.5, 7.75, -16.5, -8.0, ni=-300.0),
               Coil("mid", 4.5, 7.75, -6.5, 6.5, ni=600.0),
               Coil("hi", 4.5, 7.75, 8.0, 16.5, ni=-300.0)],
        magnets=[MagnetPart(4.0, 12.0, -6.1, magdir=90.0),   # N up
                 MagnetPart(4.0, 12.0, +6.1, magdir=270.0)], # N down (N faces N)
    )
    meta = {
        "name": "opp24c3",
        "worst_power_W": 0.364,          # 0.8 N sym, fill 0.6, hot Cu
        "worst_power_fill05_W": 0.44,
        "package_len_mm": 35.0,          # coils +/-16.5 + 2 mm end caps (limit 35)
        "base_ni": [-300.0, 600.0, -300.0],   # lo / mid / hi A-turns at s=1
        "gain_table": {                  # N per unit ni_scale, fine-sleeve run
            -2.5: 5.9697, -1.5: 6.5099, -0.5: 6.7680,
            0.5: 6.7703, 1.5: 6.5081, 2.5: 5.9711},
        "gain_source_csv": "results/opp24c3f_acd4.csv",
        "cogging_max_N": 0.0,            # true 0 (all-air); noise floor <=0.051
        "tau_us": 122.0,
        "km_worst": 1.33,                # N/sqrt(W)
        "notes": "ABC 120 vs 60: <=0.3%; fine sleeve vs automesh: <0.2%. "
                 "Best pure single-magnet dual fallback: a20 at 0.697 W.",
    }
    return model, meta


def w_o3_m12_r():
    """Finalist 2 — opposed pair + center steel washer (sim/explore_wildcard.py
    candidate "w_o3_m12_r").

    Two N52 D8x12, N poles facing a 2 mm x 8 mm-OD 1018 washer (mover);
    ironless 3-coil stator at r 4.75-7.75, current-ratio-tuned outers
    (|NI_outer/NI_mid| = 115/300 ~ power-optimal 0.38). True cogging is
    identically zero (no stator steel); washer eddy lag bounded analytically
    at 0.15-0.3 ms (FEMM harmonic solves hang headless — see the TOOLCHAIN
    FINDING in sim/explore_wildcard.py).
    """
    model = MoverModel(
        name="w_o3_m12_r",
        magnet=Magnet(radius=4.0, length=26.0, z_center=0.0),  # placeholder, unused
        coils=[Coil("mid", 4.75, 7.75, -5.0, 5.0, ni=+300.0),
               Coil("hi", 4.75, 7.75, 5.5, 15.5, ni=-115.0),
               Coil("lo", 4.75, 7.75, -15.5, -5.5, ni=-115.0)],
        magnets=[MagnetPart(4.0, 12.0, -7.0, magdir=90.0),
                 MagnetPart(4.0, 12.0, +7.0, magdir=270.0)],
        mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.0, 1.0)],
    )
    meta = {
        "name": "w_o3_m12_r",
        "worst_power_W": 0.423,          # 0.8 N sym, fill 0.6 (0.476 W from
                                         # analyze.evaluate incl. noise floor)
        "worst_power_fill05_W": 0.51,
        "package_len_mm": 33.0,          # mover 26 + stroke 5 + 2 mm caps
        "base_ni": [300.0, -115.0, -115.0],   # mid / hi / lo A-turns at s=1
        "gain_table": {                  # N per unit ni_scale, fine run
            -2.5: 2.8430, -2.0: 3.0103, -1.5: 3.1382, -1.0: 3.2309,
            -0.5: 3.2837, 0.0: 3.3014, 0.5: 3.2831, 1.0: 3.2299,
            1.5: 3.1388, 2.0: 3.0098, 2.5: 2.8409},
        "gain_source_csv": "results/w_o3_m12_r_wcgF.csv",
        "cogging_max_N": 0.0,            # true 0 (no stator steel); floor <=0.054
        "tau_us": 92.0,
        "km_worst": 1.23,                # N/sqrt(W)
        "notes": "Washer diffusion tau 0.15-0.3 ms (bias-sat mur 50-100), "
                 "3 ms at unbiased mur 1000 -> >30x slew margin.",
    }
    return model, meta


def st_w03():
    """Finalist 3 — thin steel shell (sim/explore_steel_shell.py variant
    "st_w03", refined mesh settings of its FINAL "st_w03_fine").

    Single N52 D8x12 mover, dual push-pull coil, 0.3 mm-wall 1018 tube
    r 7.7-8.0 doubling as the housing. Shortest and cheapest package; pays
    ~2.4x the hold power of opp24c3 and a small real cogging (~0.07 N,
    firmware-compensable). Pot-core end washers were REJECTED (0.4-0.9 N
    cogging). Eddy lag through the 0.3 mm wall is negligible (skin depth >
    wall to ~kHz; measured via complex flux linkage in sim/explore_steel_ac.py).

    Builder ships the FINAL fine-mesh settings (mesh_air 0.4 + 0.25 mm air
    annulus in the bore gap) — needed to resolve the ~0.07 N cogging above
    the automesh noise floor. Set mesh_air=0, fine_annulus=None to reproduce
    the faster batch-phase automesh model.
    """
    model = ShellModel(
        name="st_w03",
        magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
        coils=[Coil("lo", 4.75, 7.65, -12.5, -0.5, ni=-300.0),
               Coil("hi", 4.75, 7.65, 0.5, 12.5, ni=+300.0)],
        steels=[Steel("shell", 7.7, 8.0, -13.0, 13.0)],
        mesh_air=0.4,
        fine_annulus=(4.05, 4.7, -9.7, 9.7, 0.25),
    )
    meta = {
        "name": "st_w03",
        "worst_power_W": 1.013,          # analyze.evaluate worst direction incl.
                                         # cogging + F-I curvature, fill 0.6
        "worst_power_fill05_W": 1.22,
        "package_len_mm": 28.0,          # shell +/-13 + 2 mm end caps
        "base_ni": [-300.0, 300.0],      # lo / hi A-turns at s=1
        "gain_table": {                  # N per unit ni_scale, fine-mesh run
            -2.5: 2.2906, -2.0: 2.3784, -1.5: 2.4438, -1.0: 2.4895,
            -0.5: 2.5169, 0.0: 2.5262, 0.5: 2.5176, 1.0: 2.4896,
            1.5: 2.4435, 2.0: 2.3775, 2.5: 2.2892},
        "gain_source_csv": "results/st_w03_fine_stloptfin.csv",
        "cogging_max_N": 0.068,          # at z=+2.5, fine mesh (real, not noise)
        "tau_us": 162.0,
        "km_worst": 0.79,                # N/sqrt(W)
        "notes": "P_sym (gain-only) is 0.86 W; 1.01 W is the worse direction "
                 "with cogging. Saturation: |1 - g2/g05| < 2% at s in [-2,2].",
    }
    return model, meta


DESIGNS = {
    "opp24c3": opp24c3,
    "w_o3_m12_r": w_o3_m12_r,
    "st_w03": st_w03,
}


# =====================================================================
# CLI: re-run analyze.summarize on a design and compare force gains
# against the explore-phase reference table.
# =====================================================================

def main(argv: list[str]):
    if not argv or argv[0] not in DESIGNS:
        raise SystemExit(f"usage: designs.py <{'|'.join(DESIGNS)}> "
                         "[z1,z2,...] [s1,s2,...] [timeout_s]")
    model, meta = DESIGNS[argv[0]]()
    zs = ([float(x) for x in argv[1].split(",")] if len(argv) > 1
          else sorted(meta["gain_table"]))
    scales = ([float(x) for x in argv[2].split(",")] if len(argv) > 2
              else [-1.0, 0.0, 1.0])
    timeout = int(argv[3]) if len(argv) > 3 else 1800

    rows = run_sweep(model, zs, scales, tag="designs", timeout=timeout)
    print(summarize(evaluate(model, rows)))

    byz: dict[float, dict[float, float]] = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
    print(f"gain check vs {meta['gain_source_csv']} "
          "(gain = dF over the swept ni_scale range):")
    for z in sorted(byz):
        d = byz[z]
        smax, smin = max(d), min(d)
        gain = (d[smax] - d[smin]) / (smax - smin)
        ref = meta["gain_table"].get(z)
        if ref:
            print(f"  z={z:+5.2f}  gain={gain:8.4f}  ref={ref:8.4f}  "
                  f"diff={100 * (gain / ref - 1):+6.2f}%")
        else:
            print(f"  z={z:+5.2f}  gain={gain:8.4f}  (no reference)")
    print(f"meta: P_worst={meta['worst_power_W']} W (fill 0.5: "
          f"{meta['worst_power_fill05_W']} W), package "
          f"{meta['package_len_mm']} mm, base NI {meta['base_ni']}")


if __name__ == "__main__":
    main(sys.argv[1:])
