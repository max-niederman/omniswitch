"""Wildcard topology exploration for the omniswitch actuator (agent: wildcard).

Topologies beyond the four in sim/candidates.py:

  w_washer_dual  (b) single D8x12 N52 magnet + 1 mm steel pole washers glued
                 to each pole face (mover), aircore_dual coil pair. No stator
                 steel -> true cogging identically 0; the washers concentrate
                 pole flux radially into the coils.

  w_opposed3    (a)+(b) hybrid: two D8x10 N52 magnets, like (N) poles facing
                 a 2 mm steel washer on the mover; three-coil stack (center
                 coil sits in the outward radial flux band at the washer,
                 outer coils in the inward return bands, opposite polarity).
                 A radial-field voice coil built from commodity axial magnets.

  w_opposed3_ow  same + 1 mm outer pole washers on the far poles.

MoverModel extends femm.Model with a multi-part mover (list of MagnetPart +
mover steels that translate with z; the base class keeps mover Steel parts at
fixed absolute z, which is wrong for washers) and implements fine_box (the
base class declares but never uses it): an air rect labelled at mesh_air
wrapped around the mover so the stress-tensor integration air is finely
meshed. The fine box omits its axis edge (open_axis) to avoid overlapping
collinear segments with magnet axis edges - identical shared edges between
butted parts are deduped python-side.

All stator parts are copper/air only -> true cogging is identically zero and
air-core tau logic applies; mover-steel eddy lag is checked with an AC
harmonic force-vs-frequency sweep (acheck mode: PM Hc set to 0, coils driven,
time-average Fz vs f) plus an analytic diffusion bound.

TOOLCHAIN FINDING (2026-07-28): FEMM 4.2 HARMONIC (freq > 0) mi_analyze()
HANGS HEADLESS UNDER WINE in every configuration tested - nonlinear steel,
linear steel, and even a coils+air-only problem (see acdbg_* results; hang is
at the first "analyze" marker, geometry builds fine). The CLAUDE.md AC
eddy-lag check is therefore impossible with this toolchain for ANY variant;
steel-bearing designs must rely on analytic diffusion bounds. For this
design: washer diffusion tau = mu0*mur_inc*sigma*t^2/pi^2 = 0.15-0.3 ms at
the bias-realistic mur_inc 50-100 (washer saturated by the two magnets),
3 ms even at unbiased mur 1000 - and the primary force path (static PM field
x coil current, tau_elec = 92 us) needs no diffusion through steel at all,
so the <10 ms slew spec has >30x margin under the most pessimistic bound.

Run (from repo root, inside nix develop):
    python sim/explore_wildcard.py run <name> <tag> [fine]
    python sim/explore_wildcard.py report <name> <tag>       # re-parse csv
    python sim/explore_wildcard.py diag <name> <tag>         # per-coil gains
    python sim/explore_wildcard.py acheck <name> <tag>       # eddy-lag AC
    python sim/explore_wildcard.py abc <name> <tag>          # ABC=120 check
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from dataclasses import dataclass, field, replace

from femm import MU0, PROJECT_ROOT, RESULTS_DIR, Coil, Magnet, Model, run_sweep, z_path
from analyze import evaluate, summarize

Z_SWEEP = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
Z_FINE = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
SCALES = [-1.0, -0.5, 0.0, 0.5, 1.0]
STROKE_HALF = 2.5
CAPS_MM = 2.0
R_IN, R_OUT = 4.75, 7.75   # shared radial build: magnet r4 + bobbin/clearance,
                           # 0.25 mm housing outside -> 8.0 mm outer radius


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
    magnets: list[MagnetPart] = field(default_factory=list)
    mover_steels: list[MoverSteel] = field(default_factory=list)

    # ------------------------------------------------------------------ geometry
    # ac=True: PM Hc -> 0 and mover steel swapped for a LINEAR conductive
    # steel (mu_r = ac_steel_mu, sigma = 5.8 MS/m). Nonlinear "1018 Steel" in
    # a harmonic problem raises a modal dialog headless (verified: hang at
    # first mi_analyze), and linear is more representative anyway - the
    # washer is bias-saturated by the magnets, incremental mu_r ~ 100.
    ac_steel_mu: float = 100.0

    def _emit_doc(self, w, z: float, freq: float = 0.0, ac: bool = False):
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
        w(f'mi_probdef({freq:.6g}, "millimeters", "axi", 1e-8, 0, 30)')
        w('mi_getmaterial("Air")')
        w('mi_getmaterial("Copper")')
        for i, mp in enumerate(self.magnets):
            hc = 0.0 if ac else mp.hc
            w(f'mi_addmaterial("Magnet{i}", {mp.mu_r}, {mp.mu_r}, {hc:.1f}, '
              "0, 0.667, 0, 0, 1, 0, 0, 0)")
        if ac:
            w(f'mi_addmaterial("SteelAC", {self.ac_steel_mu:.6g}, '
              f"{self.ac_steel_mu:.6g}, 0, 0, 5.8, 0, 0, 1, 0, 0, 0)")
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
                  "SteelAC" if ac else st.material, group=1,
                  mesh=self.mesh_air)

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

    # ------------------------------------------------------------------ scripts
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

    def ac_lua(self, freqs, out_csv, fem_file, ni_scale=1.0, z=0.0) -> str:
        """Harmonic force-vs-frequency check: PM Hc -> 0, coils driven at
        base NI * ni_scale (peak), time-average Fz on the mover vs f. A flat
        |Fz(f)|/|Fz(1 Hz)| up to a few hundred Hz means mover-steel eddies do
        not add a force lag anywhere near the 10 ms slew budget."""
        L: list[str] = []
        w = L.append
        w(f'LOG = "{z_path(out_csv)}.log"')
        w("function mark(s)")
        w('    local h = openfile(LOG, "a")')
        w('    write(h, s, "\\n")')
        w("    closefile(h)")
        w("end")
        w(f'handle = openfile("{z_path(out_csv)}", "w")')
        w('write(handle, "freq,Fz\\n")')
        for f_ in freqs:
            self._emit_doc(w, z, freq=f_, ac=True)
            for c in self.coils:
                w(f'mi_modifycircprop("{c.name}", 1, {c.ni * ni_scale:.6g})')
            w(f'mark("f={f_:.6g} analyze")')
            w("mi_analyze()")
            w('mark("solved")')
            w("mi_loadsolution()")
            w("mo_groupselectblock(1)")
            w("fz = mo_blockintegral(19)")
            w(f'write(handle, {f_:.6g}, ",", fz, "\\n")')
            w("mo_close()")
            w("mi_close()")
        w("closefile(handle)")
        w("quit()")
        return "\n".join(L) + "\n"


# ---------------------------------------------------------------------- helpers
def mover_extent(model: MoverModel):
    los, his = [], []
    for mp in model.magnets:
        los.append(mp.z_offset - mp.length / 2)
        his.append(mp.z_offset + mp.length / 2)
    for st in model.mover_steels:
        los.append(st.z_bot)
        his.append(st.z_top)
    return min(los), max(his)


def package_length(model: MoverModel) -> float:
    zlo = min(c.z_bot for c in model.coils)
    zhi = max(c.z_top for c in model.coils)
    for s in model.steels:
        if not s.mover:
            zlo, zhi = min(zlo, s.z_bot), max(zhi, s.z_top)
    mlo, mhi = mover_extent(model)
    zlo = min(zlo, mlo - STROKE_HALF)
    zhi = max(zhi, mhi + STROKE_HALF)
    return (zhi - zlo) + CAPS_MM


def auto_fine_box(model: MoverModel):
    mlo, mhi = mover_extent(model)
    return (4.4, mlo - STROKE_HALF - 2.0, mhi + STROKE_HALF + 2.0)


def report(model: MoverModel, rows):
    ev = evaluate(model, rows)
    print(summarize(ev))
    byz = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
    print("per-z: gain=(F(+1)-F(-1))/2 [N/unit s], F(+1), F(-1), F(0), "
          "curv=F(+1)+F(-1)-2F(0)")
    for z in sorted(byz):
        d = byz[z]
        if 1.0 in d and -1.0 in d:
            gain = (d[1.0] - d[-1.0]) / 2
            f0 = d.get(0.0, float("nan"))
            curv = d[1.0] + d[-1.0] - 2 * d.get(0.0, 0.0)
            print(f"  z={z:+5.2f}  gain={gain:8.4f}  F+={d[1.0]:8.4f}  "
                  f"F-={d[-1.0]:8.4f}  F0={f0:8.4f}  curv={curv:+8.4f}")
    print(f"package length = {package_length(model):.1f} mm "
          f"(mover extent {mover_extent(model)}, +/-{STROKE_HALF} stroke, "
          f"+{CAPS_MM} caps)  [limit 35]")
    p1 = sum(c.power(c.ni) for c in model.coils)
    print(f"total copper power at s=1: {p1:.2f} W")
    return ev


def run_ac(model: MoverModel, freqs, tag, ni_scale=1.0, z=0.0, timeout=1200):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, f"{model.name}_{tag}")
    lua_file, out_csv, fem_file = base + ".lua", base + ".csv", base + ".fem"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua_file, "w") as f:
        f.write(model.ac_lua(freqs, out_csv, fem_file, ni_scale, z))
    subprocess.run(["femm-lua", lua_file], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    out = []
    with open(out_csv) as f:
        rdr = csv.reader(f)
        next(rdr)
        for row in rdr:
            # harmonic results are nominally real; guard against "a+I*b" text
            val = row[1].split("+I*")[0].split("+j")[0]
            out.append((float(row[0]), float(val)))
    return out


# ------------------------------------------------------------------- candidates
CANDS: dict[str, MoverModel] = {}


def _add(m: MoverModel) -> MoverModel:
    CANDS[m.name] = m
    return m


def _c(name, z1, z2, ni):
    return Coil(name, R_IN, R_OUT, z1, z2, ni=ni)


# (b) pole washers on the aircore_dual layout. D8x12 magnet, 1 mm 1018 discs
# glued to each pole. Mover 14 mm; coils as aircore_dual for an apples-to-
# apples washer-benefit readout vs the baseline (gain 2.03-2.24 N/unit).
_add(MoverModel(
    name="w_washer_dual",
    magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),   # placeholder, unused
    coils=[_c("lo", -12.5, -0.5, -300.0),
           _c("hi", 0.5, 12.5, +300.0)],
    magnets=[MagnetPart(4.0, 12.0, 0.0, magdir=90.0)],
    mover_steels=[MoverSteel("wb", 0.0, 4.0, -7.0, -6.0),
                  MoverSteel("wt", 0.0, 4.0, 6.0, 7.0)],
))

# (a)+(b) opposed pair: two D8x10, N poles facing a 2 mm center washer.
# Radial band shoots outward at the washer (+r), returns inward near the far
# poles (+/-11). Center coil + polarity, outer coils - polarity.
# Mover 22 mm -> extent 27 with stroke; stator 31 -> package 33 mm.
_add(MoverModel(
    name="w_opposed3",
    magnet=Magnet(radius=4.0, length=22.0, z_center=0.0),   # placeholder, unused
    coils=[_c("mid", -5.0, 5.0, +300.0),
           _c("hi", 5.5, 15.5, -300.0),
           _c("lo", -15.5, -5.5, -300.0)],
    magnets=[MagnetPart(4.0, 10.0, -6.0, magdir=90.0),    # N up, N pole at z=-1
             MagnetPart(4.0, 10.0, +6.0, magdir=270.0)],  # N down, N pole at z=+1
    mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.0, 1.0)],
))

# variant: + 1 mm outer pole washers (sharper return bands at the far poles)
_add(MoverModel(
    name="w_opposed3_ow",
    magnet=Magnet(radius=4.0, length=24.0, z_center=0.0),   # placeholder, unused
    coils=[_c("mid", -5.0, 5.0, +300.0),
           _c("hi", 5.5, 15.5, -300.0),
           _c("lo", -15.5, -5.5, -300.0)],
    magnets=[MagnetPart(4.0, 10.0, -6.0, magdir=90.0),
             MagnetPart(4.0, 10.0, +6.0, magdir=270.0)],
    mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.0, 1.0),
                  MoverSteel("owb", 0.0, 4.0, -12.0, -11.0),
                  MoverSteel("owt", 0.0, 4.0, 11.0, 12.0)],
))


# Iteration 2 (from diag on w_opposed3: per-coil gains mid 2.48/2.16, outer
# 0.99/1.10..0.47 N per unit at equal NI=300 and equal 4.12 W -> power-optimal
# outer/mid NI ratio r = g_out_sum/(2*g_mid) ~ 0.38 -> outers at -115 A-t).
_add(MoverModel(
    name="w_o3_r",     # ratio-tuned currents, geometry unchanged
    magnet=Magnet(radius=4.0, length=22.0, z_center=0.0),   # placeholder
    coils=[_c("mid", -5.0, 5.0, +300.0),
           _c("hi", 5.5, 15.5, -115.0),
           _c("lo", -15.5, -5.5, -115.0)],
    magnets=[MagnetPart(4.0, 10.0, -6.0, magdir=90.0),
             MagnetPart(4.0, 10.0, +6.0, magdir=270.0)],
    mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.0, 1.0)],
))

_add(MoverModel(
    name="w_o3_t3_r",  # 3 mm center washer (less saturation), ratio-tuned
    magnet=Magnet(radius=4.0, length=23.0, z_center=0.0),   # placeholder
    coils=[_c("mid", -5.5, 5.5, +300.0),
           _c("hi", 6.0, 15.5, -115.0),
           _c("lo", -15.5, -6.0, -115.0)],
    magnets=[MagnetPart(4.0, 10.0, -6.5, magdir=90.0),
             MagnetPart(4.0, 10.0, +6.5, magdir=270.0)],
    mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.5, 1.5)],
))

_add(MoverModel(
    name="w_o3_m12_r",  # D8x12 magnets (more flux), 2 mm washer, ratio-tuned
    magnet=Magnet(radius=4.0, length=26.0, z_center=0.0),   # placeholder
    coils=[_c("mid", -5.0, 5.0, +300.0),
           _c("hi", 5.5, 15.5, -115.0),
           _c("lo", -15.5, -5.5, -115.0)],
    magnets=[MagnetPart(4.0, 12.0, -7.0, magdir=90.0),
             MagnetPart(4.0, 12.0, +7.0, magdir=270.0)],
    mover_steels=[MoverSteel("cw", 0.0, 4.0, -1.0, 1.0)],
))


# --------------------------------------------------------------------------- cli
def main(argv):
    mode = argv[0]
    if mode == "run":
        name, tag = argv[1], argv[2]
        model = CANDS[name]
        fine = "fine" in argv[3:]
        if fine:
            model = replace(model, mesh_air=0.5, fine_box=auto_fine_box(model))
        zs = Z_FINE if fine else Z_SWEEP
        rows = run_sweep(model, zs, SCALES, tag=tag)
        report(model, rows)
    elif mode == "report":
        name, tag = argv[1], argv[2]
        model = CANDS[name]
        path = os.path.join(RESULTS_DIR, f"{name}_{tag}.csv")
        with open(path) as f:
            rows = [{k: float(v) for k, v in r.items()}
                    for r in csv.DictReader(f)]
        report(model, rows)
    elif mode == "diag":     # per-coil force gain (others' copper is inert)
        name, tag = argv[1], argv[2]
        model = CANDS[name]
        for c in model.coils:
            m1 = replace(model, name=f"{model.name}_only{c.name}", coils=[c])
            rows = run_sweep(m1, [-2.5, 0.0, 2.5], [-1.0, 1.0], tag=tag)
            byz = {}
            for r in rows:
                byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
            for z in sorted(byz):
                d = byz[z]
                g = (d[1.0] - d[-1.0]) / 2
                print(f"coil {c.name:4s} z={z:+5.2f} gain={g:+8.4f} N/unit "
                      f"(P(s=1)={c.power(c.ni):.2f} W)")
    elif mode == "acheck":
        name, tag = argv[1], argv[2]
        mu = float(argv[3]) if len(argv) > 3 else 100.0
        model = replace(CANDS[name], ac_steel_mu=mu)
        pts = run_ac(model, [1.0, 30.0, 100.0, 300.0, 1000.0], tag=tag)
        f0 = pts[0][1]
        print(f"steel mu_r (linear) = {mu}")
        for f_, fz in pts:
            r = fz / f0 if f0 else float("nan")
            print(f"f={f_:8.1f} Hz  Fz_avg={fz:+.6f} N  ratio_vs_1Hz={r:+.3f}")
    elif mode == "abc":
        name, tag = argv[1], argv[2]
        model = replace(CANDS[name], abc_radius=120.0)
        rows = run_sweep(model, [-2.5, 0.0], [-1.0, 1.0], tag=tag)
        byz = {}
        for r in rows:
            byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
        for z in sorted(byz):
            d = byz[z]
            print(f"ABC=120: z={z:+5.2f} gain={(d[1.0]-d[-1.0])/2:+8.4f} N/unit")
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1:])
