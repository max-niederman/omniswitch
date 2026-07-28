"""Air-core dual-coil push-pull family optimization (agent tag prefix: acd*).

Starts from `aircore_dual` in sim/candidates.py. Free parameters explored:
magnet length (D8 commodity: 8/10/12/15/20, stacked pairs incl. opposed),
coil axial split/length/gap, coil r_in (floor 4.5), total stack length.

Objective: minimize worst-case ohmic power for 0.8 N in the worse direction
anywhere in the 5 mm stroke (z_center in [-2.5, +2.5]); secondary: flatness.

Methodology notes
-----------------
* Air core => F linear in NI (verified again in round 1 with scales 0/0.5/1).
  Sweeps use ni_scales [0, 1]; the force *gain* per unit scale,
  gain(z) = F(1) - F(0), exactly cancels the constant-per-mesh weighted
  stress-tensor offset (the dominant FEMM error, ~+/-0.06 N, see the dipole
  verification), so ranking uses the noise-cancelled symmetric power
  P_sym(z) = total_power(0.8 / gain(z)).  True cogging is identically 0
  (no steel anywhere); reported s=0 rows are the numerical noise floor.
* This file has its own Lua generator (gen_lua): a superset of
  femm.Model.lua adding (a) magnet stacks [(z_offset, length, sign), ...]
  for opposed pairs, (b) an optional fine air sleeve in the magnet-coil gap
  (Model.fine_box in shared code is a no-op), (c) mover steels that actually
  translate with z.  Shared files are untouched.

Usage (inside `nix develop`, from repo root, ALWAYS under `timeout`):
    python sim/explore_aircore_dual.py list
    python sim/explore_aircore_dual.py run  <tag> <scales> <variant...>
    python sim/explore_aircore_dual.py eval <tag> <variant...>   # re-parse CSV
e.g.
    timeout 1200 nix develop . -c python sim/explore_aircore_dual.py \
        run acd1 0,1 a12 a15
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from dataclasses import dataclass, field, replace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import Coil, Magnet, Model, z_path, RESULTS_DIR, PROJECT_ROOT
from analyze import evaluate, summarize, total_power, F_TARGET


Z_SWEEP = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]


@dataclass
class Variant:
    model: Model
    # magnet stack: list of (z_offset_from_mover_center, length, sign);
    # None => single model.magnet. sign -1 => magnetized -z (magdir 270).
    stack: list | None = None
    fine_sleeve: bool = False  # 0.2 mm-meshed air sleeve in magnet-coil gap
    note: str = ""

    @property
    def mover_half_len(self) -> float:
        if self.stack is None:
            return self.model.magnet.length / 2
        lo = min(off - ln / 2 for off, ln, _ in self.stack)
        hi = max(off + ln / 2 for off, ln, _ in self.stack)
        return max(abs(lo), abs(hi))

    def package_len(self) -> float:
        zlo = min(c.z_bot for c in self.model.coils)
        zhi = max(c.z_top for c in self.model.coils)
        for s in self.model.steels:
            if not s.mover:
                zlo, zhi = min(zlo, s.z_bot), max(zhi, s.z_top)
        half = self.mover_half_len + 2.5
        return max(zhi, half) - min(zlo, -half) + 2.0  # +2 mm end caps


VARIANTS: dict[str, Variant] = {}


def _add(name: str, v: Variant) -> Variant:
    v.model = replace(v.model, name=name)  # copy: variants may share a Model
    VARIANTS[name] = v
    return v


def dual(lm: float, gap_half: float, coil_top: float, r_in: float = 4.5,
         r_out: float = 7.75, ni: float = 300.0) -> Variant:
    """Symmetric dual-coil push-pull around a single centered magnet."""
    return Variant(model=Model(
        name="tmp",
        magnet=Magnet(radius=4.0, length=lm, z_center=0.0),
        coils=[Coil("lo", r_in, r_out, -coil_top, -gap_half, ni=-ni),
               Coil("hi", r_in, r_out, gap_half, coil_top, ni=+ni)]))


# ---------------------------------------------------------------- round 1
# All r_in 4.5 (0.5 mm guide wall floor), r_out 7.75 (+0.25 housing = OD 16).
# Coil windows sized to cover each pole's +/-2.5 mm travel plus a few mm of
# radial-flux lobe on both sides.
_add("a10", dual(10, 0.5, 11.0))    # poles +/-5, travel [2.5,7.5]
_add("a12", dual(12, 0.5, 12.5))    # baseline geometry, r_in 4.75->4.5
_add("a15", dual(15, 2.0, 15.0))    # poles +/-7.5, travel [5,10]
_add("a20", dual(20, 4.5, 16.0))    # poles +/-10, travel [7.5,12.5]; pkg 34

# ---------------------------------------------------------------- round 2
# Round-1 result: P_worst(sym) a20 0.697 < a15 0.785 < a12 0.922 < a10 1.089 W
# -> longer magnet wins; tune a20 gap/length + test opposed-pair topology.
_add("a20g3", dual(20, 3.0, 16.5))   # more inner copper, pkg 35.0 (limit)
_add("a20s", dual(20, 5.5, 16.0))    # trimmed inner copper, pkg 34
_add("a15g1", dual(15, 1.0, 15.0))   # smaller mid-gap, +1 mm copper each
_add("a15g3", dual(15, 3.0, 15.0))   # bigger mid-gap (skip B_r sign-flip zone)
_add("a15s", dual(15, 2.0, 13.5))    # trimmed outer copper (dilution check)
_add("a15L", dual(15, 2.0, 16.5))    # stretched outer copper, pkg 35.0
# Opposed pair: 2x D8x10 N-to-N (0.2 mm spacer), single fat center coil
# capturing the doubled radial flux at the joint.
_add("opp20", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=10.0, z_center=0.0),
                coils=[Coil("main", 4.5, 7.75, -6.0, 6.0, ni=600.0)]),
    stack=[(-5.1, 10.0, +1), (+5.1, 10.0, -1)],
    note="opposed pair D8x10+D8x10, center coil"))

# ---------------------------------------------------------------- round 3
# Round-2 result: opp20c3 0.477 W << a20* ~0.70 W. Reciprocity check
# (F_i = NI_i dlam_i/dz on the s=0 flux columns) reproduces the total gain
# and shows the -300/600/-300 split is within 1-3% of the optimal current
# split -> tune geometry, not ratio: extend outer coils into the package
# headroom (31 -> 35 mm) and try longer opposed stacks (2x D8x12; 2x D8x15
# would need pkg 37.2 mm -> does not fit).
_add("opp20c3L", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=10.0, z_center=0.0),
                coils=[Coil("lo", 4.5, 7.75, -16.5, -6.5, ni=-300.0),
                       Coil("mid", 4.5, 7.75, -6.0, 6.0, ni=600.0),
                       Coil("hi", 4.5, 7.75, 6.5, 16.5, ni=-300.0)]),
    stack=[(-5.1, 10.0, +1), (+5.1, 10.0, -1)],
    note="opp pair D8x10, outer coils extended to +/-16.5, pkg 35"))
_add("opp24c3", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
                coils=[Coil("lo", 4.5, 7.75, -16.5, -8.0, ni=-300.0),
                       Coil("mid", 4.5, 7.75, -6.5, 6.5, ni=600.0),
                       Coil("hi", 4.5, 7.75, 8.0, 16.5, ni=-300.0)]),
    stack=[(-6.1, 12.0, +1), (+6.1, 12.0, -1)],
    note="opp pair D8x12, outer poles +/-12.2, pkg 35"))
_add("opp20w3", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=10.0, z_center=0.0),
                coils=[Coil("lo", 4.5, 7.75, -16.5, -8.0, ni=-300.0),
                       Coil("mid", 4.5, 7.75, -7.5, 7.5, ni=600.0),
                       Coil("hi", 4.5, 7.75, 8.0, 16.5, ni=-300.0)]),
    stack=[(-5.1, 10.0, +1), (+5.1, 10.0, -1)],
    note="opp pair D8x10, wide mid coil +/-7.5, pkg 35"))
_add("opp24w3", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=12.0, z_center=0.0),
                coils=[Coil("lo", 4.5, 7.75, -16.5, -8.0, ni=-300.0),
                       Coil("mid", 4.5, 7.75, -7.5, 7.5, ni=600.0),
                       Coil("hi", 4.5, 7.75, 8.0, 16.5, ni=-300.0)]),
    stack=[(-6.1, 12.0, +1), (+6.1, 12.0, -1)],
    note="opp pair D8x12, wide mid coil +/-7.5, pkg 35"))
_add("opp20c3", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=10.0, z_center=0.0),
                coils=[Coil("lo", 4.5, 7.75, -14.5, -6.5, ni=-300.0),
                       Coil("mid", 4.5, 7.75, -6.0, 6.0, ni=600.0),
                       Coil("hi", 4.5, 7.75, 6.5, 14.5, ni=-300.0)]),
    stack=[(-5.1, 10.0, +1), (+5.1, 10.0, -1)],
    note="opposed pair + outer bucking coils, pkg 31+2"))
_add("opp16", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=8.0, z_center=0.0),
                coils=[Coil("main", 4.5, 7.75, -6.0, 6.0, ni=600.0)]),
    stack=[(-4.1, 8.0, +1), (+4.1, 8.0, -1)],
    note="opposed pair D8x8, center coil"))
_add("opp20w", Variant(
    model=Model(name="tmp",
                magnet=Magnet(radius=4.0, length=10.0, z_center=0.0),
                coils=[Coil("main", 4.5, 7.75, -7.5, 7.5, ni=600.0)]),
    stack=[(-5.1, 10.0, +1), (+5.1, 10.0, -1)],
    note="opposed pair D8x10, wider center coil"))
_add("a15g1s", dual(15, 1.0, 14.0))  # small gap + mild outer trim

# ---------------------------------------------------------------- round 4
# Round-3 result: opp24c3 0.364 W < opp24w3 0.398 < opp20c3L 0.455 <
# opp20c3 0.477 < opp20w3 0.573 W.  Winner = opp24c3.  Validation round:
# fine-sleeve mesh + bipolar scales, ABC-radius doubling, denser z-grid.
_add("opp20c3f", replace(VARIANTS["opp20c3"],
                         model=VARIANTS["opp20c3"].model,
                         fine_sleeve=True))
_add("a15g1f", replace(VARIANTS["a15g1"], fine_sleeve=True))
_add("opp24c3f", replace(VARIANTS["opp24c3"], fine_sleeve=True))
_add("opp24c3r120", replace(
    VARIANTS["opp24c3"],
    model=replace(VARIANTS["opp24c3"].model, abc_radius=120.0)))

# Round-4 validation of opp24c3 (tags acd4/acd4z, results/ opp24c3*_acd4*):
# * fine sleeve + scales -1/0/1: gains 5.970/6.510/6.768/6.770/6.508/5.971
#   vs automesh 5.967/6.516/6.779/6.778/6.523/5.980 -> <0.2% shift; noise
#   floor F0 drops to <=0.051 N; s_up/s_dn symmetric after F0 subtraction.
# * ABC 120 vs 60 mm: gains match to <=0.3% -> boundary converged.
# * dense z (-2,0,2): 6.283/6.807/6.278 -> smooth profile, no dips between
#   the 6 sampled stroke positions.
# FINAL (opp24c3): P_worst(0.8 N, worse direction, fill 0.6, hot Cu)
#   = 0.364 W (0.44 W if fill derated to 0.5); Km_worst = 1.33 N/sqrt(W);
#   package 35.0 mm; tau 122 us; true cogging 0 (all-air).  Best pure
#   single-magnet dual fallback: a20 (D8x20) at 0.697 W, package 34 mm.


# ---------------------------------------------------------------- FEMM I/O

def gen_lua(v: Variant, z_positions: list[float], ni_scales: list[float],
            out_csv: str, fem_file: str) -> str:
    """Superset of femm.Model.lua: magnet stacks, fine sleeve, moving steel."""
    model = v.model
    m = model.magnet
    stack = v.stack or [(0.0, m.length, +1)]
    L: list[str] = []
    w = L.append

    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')
    flux_cols = "".join(f",flux_{c.name}" for c in model.coils)
    w(f'write(handle, "z,ni_scale,Fz{flux_cols}\\n")')

    def rect(r1, z1, r2, z2):
        for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                             (r2, z2, r1, z2), (r1, z2, r1, z1)]:
            w(f"mi_addnode({a:.6g}, {b:.6g})")
            w(f"mi_addnode({c:.6g}, {d:.6g})")
            w(f"mi_addsegment({a:.6g}, {b:.6g}, {c:.6g}, {d:.6g})")

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
        for s in model.steels:
            w(f'mi_getmaterial("{s.material}")')

        for (off, ln, sgn) in stack:
            zc = z + off
            rect(0, zc - ln / 2, m.radius, zc + ln / 2)
            label(m.radius / 2, zc, "Magnet",
                  magdir=(90 if sgn > 0 else 270), group=1,
                  mesh=model.mesh_air)

        for c in model.coils:
            rect(c.r_in, c.z_bot, c.r_out, c.z_top)
            w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
            label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
                  circuit=c.name, turns=1, mesh=model.mesh_air)

        for s in model.steels:
            dz = z if s.mover else 0.0
            rect(s.r_in, s.z_bot + dz, s.r_out, s.z_top + dz)
            label((s.r_in + s.r_out) / 2, (s.z_bot + s.z_top) / 2 + dz,
                  s.material, group=1 if s.mover else 0, mesh=model.mesh_air)

        if v.fine_sleeve:
            half = v.mover_half_len
            rin_min = min(c.r_in for c in model.coils)
            r1, r2 = m.radius + 0.05, rin_min - 0.05
            rect(r1, z - half - 4.0, r2, z + half + 4.0)
            label((r1 + r2) / 2, z + half + 3.5, "Air", mesh=0.2)

        w(f"mi_makeABC(7, {model.abc_radius:.6g}, 0, 0, 0)")
        label(model.abc_radius * 0.5, model.abc_radius * 0.5, "Air")

        w(f'mark("z={z:.6g} built")')
        w(f'mi_saveas("{z_path(fem_file)}")')
        for s in ni_scales:
            for c in model.coils:
                w(f'mi_modifycircprop("{c.name}", 1, {c.ni * s:.6g})')
            w(f'mark("z={z:.6g} s={s:.6g} analyze")')
            w("mi_analyze()")
            w('mark("solved")')
            w("mi_loadsolution()")
            w("mo_groupselectblock(1)")
            w("fz = mo_blockintegral(19)")
            w(f'write(handle, {z:.6g}, ",", {s:.6g}, ",", fz)')
            for c in model.coils:
                w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
                w('write(handle, ",", lam)')
            w('write(handle, "\\n")')
            w("mo_close()")
        w("mi_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def run(v: Variant, z_positions, ni_scales, tag: str, timeout: int = 1700):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, f"{v.model.name}_{tag}")
    lua_file, out_csv, fem_file = base + ".lua", base + ".csv", base + ".fem"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua_file, "w") as f:
        f.write(gen_lua(v, z_positions, ni_scales, out_csv, fem_file))
    subprocess.run(["femm-lua", lua_file], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    with open(out_csv) as f:
        return [{k: float(val) for k, val in row.items()}
                for row in csv.DictReader(f)]


# ---------------------------------------------------------------- reporting

def report(v: Variant, rows):
    print(summarize(evaluate(v.model, rows)))
    byz: dict[float, dict[float, float]] = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
    zs = sorted(byz)
    print("  noise-cancelled:  z   gain[N/u]  F0[N]   s_sym   P_sym[W]")
    worst, gains = 0.0, []
    for z in zs:
        d = byz[z]
        smax, smin = max(d), min(d)
        gain = (d[smax] - d[smin]) / (smax - smin)
        f0 = d.get(0.0, float("nan"))
        s_sym = F_TARGET / gain
        p = total_power(v.model, s_sym)
        worst = max(worst, p)
        gains.append(gain)
        print(f"    {z:6.2f}  {gain:8.3f}  {f0:7.3f}  {s_sym:6.3f}  {p:7.3f}")
    flat = (max(gains) - min(gains)) / max(gains) * 100
    print(f"  P_worst(sym 0.8N) = {worst:.3f} W ; "
          f"gain flatness (max-min)/max = {flat:.1f}%")
    print(f"  package length = {v.package_len():.1f} mm "
          f"(mover half-len {v.mover_half_len:.1f}) ; note: {v.note}")
    print()


def main(argv: list[str]):
    cmd = argv[0]
    if cmd == "list":
        for n, v in VARIANTS.items():
            print(n, "-", v.note or "dual push-pull",
                  f"pkg={v.package_len():.1f}mm")
        return
    if cmd == "run":
        tag, scales_s = argv[1], argv[2]
        scales = [float(x) for x in scales_s.split(",")]
        for name in argv[3:]:
            v = VARIANTS[name]
            rows = run(v, Z_SWEEP, scales, tag)
            report(v, rows)
        return
    if cmd == "runz":  # explicit z list: runz <tag> <scales> <zlist> <name...>
        tag, scales_s, z_s = argv[1], argv[2], argv[3]
        scales = [float(x) for x in scales_s.split(",")]
        zlist = [float(x) for x in z_s.split(",")]
        for name in argv[4:]:
            v = VARIANTS[name]
            rows = run(v, zlist, scales, tag)
            report(v, rows)
        return
    if cmd == "eval":
        tag = argv[1]
        for name in argv[2:]:
            v = VARIANTS[name]
            with open(os.path.join(RESULTS_DIR,
                                   f"{name}_{tag}.csv")) as f:
                rows = [{k: float(val) for k, val in row.items()}
                        for row in csv.DictReader(f)]
            report(v, rows)
        return
    raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main(sys.argv[1:])
