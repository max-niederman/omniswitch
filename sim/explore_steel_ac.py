"""Eddy-lag (magnetic diffusion) check for steel-return variants.

FEMM time-harmonic solves (mi_probdef freq > 0) of a steel-shell actuator with
the magnet region replaced by a LINEAR mur=1.05, Hc=0 material (permanent
magnets are invalid in harmonic solves), coils driven at the base NI, at
f = 0 (DC ref), 20, 50, 100, 200, 500 Hz.  Metric: complex single-turn flux
linkage of each coil circuit from mo_getcircuitproperties — |lam(f)|/|lam(0)|
and phase give the field-diffusion lag; fit an equivalent first-order pole.

Complex-value conventions (verified empirically, scratchpad ac_probe2):
  - FEMM's Lua 4.0 is complex-valued; abs(), arg(), conj(), I exist
    (Re()/Im() do NOT); write() prints complex as "a+I*b".
  - Nonlinear materials (matlib "1018 Steel") in a harmonic solve raise a
    modal dialog headless -> hang. Use LINEAR steel: mur=529 (1018 initial,
    conservative = most screening) and mur=100 (saturated-bias sensitivity).

Configs:
  ctrl_air : no steel, magnet sigma=0   -> validates method, expect ~0 lag
  mag_eddy : no steel, magnet sigma=0.667 MS/m (solid NdFeB eddies alone)
  steel529 : linear steel mur=529 sigma=5.8, magnet sigma=0.667
  steel100 : linear steel mur=100 sigma=5.8, magnet sigma=0.667
  slit529  : steel529 + shell in a parallel-connected 0 A circuit — kills the
             net circumferential (shorted-turn) eddy mode = axially SLIT tube
             (manual: zero-net-current parallel circuit enforces connectivity;
             only locally circulating eddies remain)

Run inside `nix develop`:
    python sim/explore_steel_ac.py <variant-name>   # e.g. st_pot05
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import PROJECT_ROOT, RESULTS_DIR, z_path  # noqa: E402
from explore_steel_shell import ShellModel, VARIANTS  # noqa: E402

TAG = "stlac"
FREQS = [0.0, 20.0, 50.0, 100.0, 200.0, 500.0]
CONFIGS = [
    # (name, include_steel, steel_mur, magnet_sigma, slit)
    ("ctrl_air", False, 0.0, 0.0, False),
    ("mag_eddy", False, 0.0, 0.667, False),
    ("steel529", True, 529.0, 0.667, False),
    ("steel100", True, 100.0, 0.667, False),
    ("slit529", True, 529.0, 0.667, True),
]


def ac_lua(m: ShellModel, out_csv: str, fem_file: str) -> str:
    """Harmonic sweep Lua: configs x freqs at magnet z_center = 0."""
    mag = m.magnet
    L: list[str] = []
    w = L.append

    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w("    local h = openfile(LOG, \"a\")")
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')
    cols = "".join(
        f",absI_{c.name},argI_{c.name},absV_{c.name},argV_{c.name}"
        f",absLam_{c.name},argLam_{c.name}" for c in m.coils)
    w(f'write(handle, "config,freq{cols}\\n")')

    def seg(r1, z1, r2, z2):
        w(f"mi_addnode({r1:.6g}, {z1:.6g})")
        w(f"mi_addnode({r2:.6g}, {z2:.6g})")
        w(f"mi_addsegment({r1:.6g}, {z1:.6g}, {r2:.6g}, {z2:.6g})")

    def rect(r1, z1, r2, z2):
        for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                             (r2, z2, r1, z2), (r1, z2, r1, z1)]:
            seg(a, b, c, d)

    def label(r, z, mat, circuit="", group=0, turns=0):
        w(f"mi_addblocklabel({r:.6g}, {z:.6g})")
        w(f"mi_selectlabel({r:.6g}, {z:.6g})")
        w(f'mi_setblockprop("{mat}", 1, 0, "{circuit}", 0, {group}, {turns})')
        w("mi_clearselected()")

    for (cfg, with_steel, mur, msig, slit) in CONFIGS:
        for f in FREQS:
            w("newdocument(0)")
            w(f'mi_probdef({f:.6g}, "millimeters", "axi", 1e-8, 0, 30)')
            w('mi_getmaterial("Air")')
            w('mi_addmaterial("CuAC", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)')
            w(f'mi_addmaterial("MagnetAC", {mag.mu_r}, {mag.mu_r}, 0, 0, '
              f"{msig:.6g}, 0, 0, 1, 0, 0, 0)")
            if with_steel:
                w(f'mi_addmaterial("SteelLin", {mur:.6g}, {mur:.6g}, 0, 0, '
                  "5.8, 0, 0, 1, 0, 0, 0)")
            steel_circ = ""
            if with_steel and slit:
                # parallel (circtype 0) circuit, 0 A net: slit-tube equivalent
                w('mi_addcircprop("eddy0", 0, 0)')
                steel_circ = "eddy0"

            # magnet at mid-stroke, linear, Hc=0 (mover geometry kept for the
            # eddy path it provides)
            rect(0, -mag.length / 2, mag.radius, mag.length / 2)
            label(mag.radius / 2, 0, "MagnetAC", group=1)

            for c in m.coils:
                rect(c.r_in, c.z_bot, c.r_out, c.z_top)
                w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
                label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "CuAC",
                      circuit=c.name, turns=1)

            if with_steel:
                for s in m.steels:
                    rect(s.r_in, s.z_bot, s.r_out, s.z_top)
                    label((s.r_in + s.r_out) / 2, (s.z_bot + s.z_top) / 2,
                          "SteelLin", circuit=steel_circ)
                for (name, pts, (lr, lz), mat) in m.steel_polys:
                    n = len(pts)
                    for i in range(n):
                        r1, z1 = pts[i]
                        r2, z2 = pts[(i + 1) % n]
                        seg(r1, z1, r2, z2)
                    label(lr, lz, "SteelLin", circuit=steel_circ)

            w(f"mi_makeABC(7, {m.abc_radius:.6g}, 0, 0, 0)")
            label(m.abc_radius * 0.5, m.abc_radius * 0.5, "Air")

            w(f'mark("{cfg} f={f:.6g} built")')
            w(f'mi_saveas("{z_path(fem_file)}")')
            w("mi_analyze()")
            w(f'mark("{cfg} f={f:.6g} solved")')
            w("mi_loadsolution()")
            w(f'write(handle, "{cfg}", ",", {f:.6g})')
            for c in m.coils:
                w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
                w('write(handle, ",", abs(ic), ",", arg(ic), ",", abs(vc), '
                  '",", arg(vc), ",", abs(lam), ",", arg(lam))')
            w('write(handle, "\\n")')
            w("mo_close()")
            w("mi_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def run_ac(m: ShellModel, timeout: int = 1800):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, f"{m.name}_{TAG}")
    lua_file, out_csv, fem_file = base + ".lua", base + ".csv", base + ".fem"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua_file, "w") as f:
        f.write(ac_lua(m, out_csv, fem_file))
    subprocess.run(["femm-lua", lua_file], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    with open(out_csv) as f:
        return list(csv.DictReader(f))


def analyze_ac(m: ShellModel, rows):
    print(f"== eddy-lag (harmonic) check: {m.name} ==")
    coil = m.coils[-1].name   # coils are symmetric; report one + cross-check
    by_cfg: dict[str, list] = {}
    for r in rows:
        by_cfg.setdefault(r["config"], []).append(r)
    for cfg, rr in by_cfg.items():
        rr.sort(key=lambda r: float(r["freq"]))
        ref = next(r for r in rr if float(r["freq"]) == 0.0)
        lam0 = float(ref[f"absLam_{coil}"])
        arg0 = float(ref[f"argLam_{coil}"])
        print(f"-- {cfg} (|lam(0)| = {lam0:.4e} Wb, arg0 = {arg0:.3f} rad)")
        print("   f[Hz]   |lam|/|lam0|   phase[deg]   tau_eq[ms]=tan(-phi)/w")
        for r in rr:
            f = float(r["freq"])
            rel = float(r[f"absLam_{coil}"]) / lam0
            # phase relative to DC (arg0 is 0 or pi depending on sign)
            phi = float(r[f"argLam_{coil}"]) - arg0
            while phi > math.pi:
                phi -= 2 * math.pi
            while phi < -math.pi:
                phi += 2 * math.pi
            tau = math.tan(-phi) / (2 * math.pi * f) * 1e3 if f > 0 else 0.0
            print(f"  {f:6.0f}   {rel:10.4f}   {math.degrees(phi):9.3f}   "
                  f"{tau:8.3f}")
    print()


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "st_pot05"
    m = VARIANTS[name]
    rows = run_ac(m)
    analyze_ac(m, rows)
