"""Steel-return family optimization, starting from steel_shell_dual.

Free parameters explored: shell wall thickness (0.3 / 0.5 / 1.0 mm), shell
length/overhang, end washers (pot-core style end plates), coil dimensions.
Magnet fixed at N52 D8x12 (commodity), moving-magnet armature.

Steel washers that touch the shell are modeled as a single C-shaped polygon
(one connected steel region) to avoid the coincident-collinear-segment FEMM
pitfall; coils keep a 0.05 mm air gap to all steel.

Also provides an optional fine-mesh air annulus in the bore gap (where the
weighted-stress-tensor weighting gradient lives) to knock down the ~0.05-0.1 N
automesh force-noise floor for the final cogging profile.

Run inside `nix develop`:
    python sim/explore_steel_shell.py batch    # 6-variant ranking sweep
    python sim/explore_steel_shell.py final    # refined run of the winner
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import Coil, Magnet, Model, Steel, z_path, run_sweep  # noqa: E402
from analyze import evaluate, summarize, total_power  # noqa: E402

TAG = "stlopt"          # unique agent tag (files: results/<name>_<tag>.*)
Z_SWEEP = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
Z_FINE = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
NI_SCALES = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]   # curvature-resolving


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


def pot_poly(r_bore, r_shell_in, z_wash_in, z_end):
    """C-shaped shell + two end washers as one closed polygon (r-z plane)."""
    return [
        (r_bore, -z_end), (8.0, -z_end), (8.0, z_end), (r_bore, z_end),
        (r_bore, z_wash_in), (r_shell_in, z_wash_in),
        (r_shell_in, -z_wash_in), (r_bore, -z_wash_in),
    ]


MAG = Magnet(radius=4.0, length=12.0, z_center=0.0)

VARIANTS: dict[str, ShellModel] = {}


def _add(m: ShellModel) -> ShellModel:
    VARIANTS[m.name] = m
    return m


# 0.3 mm wall, plain tube (thin: less eddy lag, more copper room)
_add(ShellModel(
    name="st_w03",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 7.65, -12.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 7.65, 0.5, 12.5, ni=+300.0)],
    steels=[Steel("shell", 7.7, 8.0, -13.0, 13.0)],
))

# 1.0 mm wall, plain tube (saturation-proof but eats copper)
_add(ShellModel(
    name="st_w10",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 6.95, -12.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 6.95, 0.5, 12.5, ni=+300.0)],
    steels=[Steel("shell", 7.0, 8.0, -13.0, 13.0)],
))

# 0.5 mm wall + 1.0 mm end washers (pot-style), same coil as baseline
_add(ShellModel(
    name="st_pot05",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 7.45, -12.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 7.45, 0.5, 12.5, ni=+300.0)],
    steel_polys=[("pot", pot_poly(4.75, 7.5, 12.55, 13.55), (7.75, 0.0),
                  "1018 Steel")],
))

# 0.3 mm wall + washers
_add(ShellModel(
    name="st_pot03",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 7.65, -12.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 7.65, 0.5, 12.5, ni=+300.0)],
    steel_polys=[("pot", pot_poly(4.75, 7.7, 12.55, 13.55), (7.85, 0.0),
                  "1018 Steel")],
))

# 0.5 mm wall, longer coils (use length budget), no washers
_add(ShellModel(
    name="st_w05L",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 7.45, -14.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 7.45, 0.5, 14.5, ni=+300.0)],
    steels=[Steel("shell", 7.5, 8.0, -15.0, 15.0)],
))

# 0.5 mm wall + washers + longer coils
_add(ShellModel(
    name="st_pot05L",
    magnet=MAG,
    coils=[Coil("lo", 4.75, 7.45, -14.5, -0.5, ni=-300.0),
           Coil("hi", 4.75, 7.45, 0.5, 14.5, ni=+300.0)],
    steel_polys=[("pot", pot_poly(4.75, 7.5, 14.55, 15.55), (7.75, 0.0),
                  "1018 Steel")],
))


def stator_extent(m: ShellModel):
    zs = []
    for c in m.coils:
        zs += [c.z_bot, c.z_top]
    for s in m.steels:
        zs += [s.z_bot, s.z_top]
    for (_, pts, _, _) in m.steel_polys:
        zs += [p[1] for p in pts]
    return min(zs), max(zs)


def package_len(m: ShellModel, stroke_half=2.5):
    lo, hi = stator_extent(m)
    mlo = -stroke_half - m.magnet.length / 2
    mhi = stroke_half + m.magnet.length / 2
    return max(hi, mhi) - min(lo, mlo) + 2.0   # +2 mm end caps


def report(m: ShellModel, rows):
    ev = evaluate(m, rows)
    print(summarize(ev))
    byz = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]
    print("  flatness/saturation (cogging-cancelling gains):")
    print("  z[mm]   g05=(F(.5)-F(-.5))/1   g1=(F(1)-F(-1))/2   "
          "g2=(F(2)-F(-2))/4   sat%=1-g2/g05")
    for z in sorted(byz):
        f = byz[z]
        g05 = (f[0.5] - f[-0.5]) / 1.0
        g1 = (f[1.0] - f[-1.0]) / 2.0
        g2 = (f[2.0] - f[-2.0]) / 4.0
        print(f"  {z:6.2f}  {g05:8.4f}  {g1:8.4f}  {g2:8.4f}  "
              f"{100 * (1 - g2 / g05):6.2f}%")
    plen = package_len(m)
    print(f"  package length: {plen:.1f} mm ({'PASS' if plen <= 35 else 'FAIL'} <= 35)")
    print()
    return ev


def main(mode: str):
    if mode == "batch":
        for name, m in VARIANTS.items():
            rows = run_sweep(m, Z_SWEEP, NI_SCALES, tag=TAG)
            report(m, rows)
    elif mode == "final":
        m = FINAL
        rows = run_sweep(m, Z_FINE, NI_SCALES, tag=TAG + "fin")
        report(m, rows)
    else:
        raise SystemExit(f"unknown mode {mode}")


# Winner (after the batch run): st_w03 — 0.3 mm plain shell, dual coil.
# Washer (pot) variants rejected: real ~0.7-0.9 N antisymmetric cogging
# (negative spring, ~100% of the 0.8 N target) at the stroke ends.
# Refined mesh for the final numbers: sized labels + fine air annulus in the
# bore gap (the WST weighting-gradient region) to cut the automesh noise floor.
import dataclasses

FINAL: ShellModel = dataclasses.replace(
    VARIANTS["st_w03"],
    name="st_w03_fine",
    mesh_air=0.4,
    fine_annulus=(4.05, 4.7, -9.7, 9.7, 0.25),
)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "batch")
