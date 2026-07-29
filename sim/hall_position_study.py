"""Hall-sensor placement study for st_w10: does board-level mounting work?

Question (2026-07-28): the crosstalk study assumed the Hall sensing element
in-bore at z = -9.5. Could it instead sit lower -- e.g. on the main PCB
~5 mm below the shell (z ~ -18) -- and still give a usable position signal?

Answer (recorded in docs/DESIGN.md "Position sensor"): board-level REJECTED
as primary (no shell shielding below the package: ~85x the in-bore neighbor
pickup; adversarial multi-key sum breaches the 2-5% force-fidelity bar);
frozen setup is the cap-riding in-bore stack, element z ~ -9.9.

Method: axisymmetric FEMM (st_w10 exact geometry from converge_shell), mover
swept over the +/-2.5 mm stroke at drives {0, +1, -1}. At each candidate
sensor height z_s we average Bz through
  * an on-axis disc r = 0..0.75  (own position signal), and
  * annuli around r = 19.05 / 26.94 / 38.10  (field at an unshielded
    neighbor's sensor, nearest/diagonal/second-ring at 19.05 mm key pitch --
    upper bound, same methodology as qualify_crosstalk)
via mo_lineintegral(0) contour windows (point evaluation hangs under Wine).

Validation per project discipline:
  valA  bare magnet (no coils/shell) vs the exact on-axis solenoid-sheet
        formula  Bz = (Br/2/mu_r)*[c(z-zA) - c(z-zB)], c(u)=u/sqrt(u^2+R^2)
  valB  ABC radius 60 -> 120 at (s=0, drive 0)
  valC  sensor-strip mesh halved at (s=0, drive 0)

Run from the repo root inside nix develop, ALWAYS under timeout:
    timeout 1700 nix develop . -c python sim/hall_position_study.py run
    python sim/hall_position_study.py analyze   # offline re-analysis
Outputs: results/hallpos.{lua,fem,csv,csv.log}
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from femm import PROJECT_ROOT, RESULTS_DIR, z_path  # noqa: E402
from converge_shell import VARIANTS  # noqa: E402  (frozen geometry source)

M = VARIANTS["st_w10"]
MAG = M.magnet                      # r=4, len=12, Br 1.43 (N52 sim convention)
SHELL = M.steels[0]                 # 7.0..8.0, z +/-13

MOVERS = [-2.5, -1.25, 0.0, 1.25, 2.5]
DRIVES = [0.0, 1.0, -1.0]           # ni_scale; +/-300 A-turns per coil
SENSOR_Z = [-9.5, -11.0, -12.0, -13.0, -14.0, -15.0, -16.5, -18.0, -20.0]
DISC_R = 0.75                       # on-axis averaging disc radius, mm
# annuli around the victim-sensor radii: nearest / diagonal / second-ring
# neighbor at 19.05 mm key pitch
NBR_WINDOWS = [("n1", 19.05, (18.55, 19.55)),
               ("n2", 26.94, (26.44, 27.44)),
               ("n3", 38.10, (37.60, 38.60))]
BASE = os.path.join(RESULTS_DIR, "hallpos")

N45SH_SCALE = 1.35 / 1.43           # linear Br rescale for reporting


def emit_doc(w, segs: set, mover_z: float, with_shell: bool, with_coils: bool,
             abc: float, mesh_scale: float = 1.0):
    def seg(r1, z1, r2, z2):
        key = tuple(sorted([(round(r1, 6), round(z1, 6)),
                            (round(r2, 6), round(z2, 6))]))
        if key in segs:
            return
        segs.add(key)
        w(f"mi_addnode({r1:.6g}, {z1:.6g})")
        w(f"mi_addnode({r2:.6g}, {z2:.6g})")
        w(f"mi_addsegment({r1:.6g}, {z1:.6g}, {r2:.6g}, {z2:.6g})")

    def rect(r1, z1, r2, z2):
        for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                             (r2, z2, r1, z2), (r1, z2, r1, z1)]:
            seg(a, b, c, d)

    def label(r, z, mat, circuit="", magdir=0, group=0, turns=0, mesh=0.0):
        automesh = 1 if mesh <= 0 else 0
        w(f"mi_addblocklabel({r:.6g}, {z:.6g})")
        w(f"mi_selectlabel({r:.6g}, {z:.6g})")
        w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, "{circuit}", '
          f"{magdir}, {group}, {turns})")
        w("mi_clearselected()")

    w("newdocument(0)")
    w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w('mi_getmaterial("Copper")')
    w(f'mi_addmaterial("Magnet", {MAG.mu_r}, {MAG.mu_r}, {MAG.hc:.1f}, '
      "0, 0.667, 0, 0, 1, 0, 0, 0)")
    w('mi_getmaterial("1018 Steel")')

    rect(0, mover_z - MAG.length / 2, MAG.radius, mover_z + MAG.length / 2)
    label(MAG.radius / 2, mover_z, "Magnet", magdir=90, group=1, mesh=0.4)

    if with_coils:
        for c in M.coils:
            rect(c.r_in, c.z_bot, c.r_out, c.z_top)
            w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
            label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
                  circuit=c.name, turns=1, mesh=0.75)

    if with_shell:
        rect(SHELL.r_in, SHELL.z_bot, SHELL.r_out, SHELL.z_top)
        label((SHELL.r_in + SHELL.r_out) / 2, 0, "1018 Steel", group=0,
              mesh=0.3)

    # fine air strips over the measurement windows (top clears the magnet
    # bottom face, which reaches z = -8.5 at s = -2.5)
    rect(0, -21.0, 2.0, -8.9)
    label(1.0, -20.5, "Air", mesh=0.3 * mesh_scale)
    rect(17.8, -21.0, 39.4, -12.0)
    label(28.0, -20.5, "Air", mesh=0.6 * mesh_scale)

    w(f"mi_makeABC(7, {abc:.6g}, 0, 0, 0)")
    label(abc * 0.5, abc * 0.5, "Air")


def emit_measure(w, tag: str, s: float, drive: float):
    for zs in SENSOR_Z:
        w("mo_clearcontour()")
        w(f"mo_addcontour(0, {zs:.6g})")
        w(f"mo_addcontour({DISC_R:.6g}, {zs:.6g})")
        w("ta, aa = mo_lineintegral(0)")
        vals = []
        for nm, _, (r1, r2) in NBR_WINDOWS:
            w("mo_clearcontour()")
            w(f"mo_addcontour({r1:.6g}, {zs:.6g})")
            w(f"mo_addcontour({r2:.6g}, {zs:.6g})")
            w(f"t{nm}, a{nm} = mo_lineintegral(0)")
            vals.append(f'",", t{nm}, ",", a{nm}')
        w("mo_clearcontour()")
        w(f'write(handle, "{tag},{s:.6g},{drive:.6g},{zs:.6g},", '
          f'ta, ",", aa, {", ".join(vals)}, "\\n")')


def gen_lua() -> str:
    L: list[str] = []
    w = L.append
    w(f'LOG = "{z_path(BASE + ".csv")}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(BASE + ".csv")}", "w")')
    ncols = "".join(f",tot_{nm},avg_{nm}" for nm, _, _ in NBR_WINDOWS)
    w(f'write(handle, "tag,s,drive,zs,tot_axis,avg_axis{ncols}\\n")')

    # main sweep: full st_w10, 5 mover positions x 3 drives
    for s in MOVERS:
        segs: set = set()
        emit_doc(w, segs, s, with_shell=True, with_coils=True, abc=60.0)
        w(f'mark("s={s:.6g} built")')
        w(f'mi_saveas("{z_path(BASE + ".fem")}")')
        for d in DRIVES:
            for c in M.coils:
                w(f'mi_modifycircprop("{c.name}", 1, {c.ni * d:.6g})')
            w(f'mark("s={s:.6g} d={d:.6g} analyze")')
            w("mi_analyze()")
            w('mark("solved")')
            w("mi_loadsolution()")
            emit_measure(w, "main", s, d)
            w("mo_close()")
        w("mi_close()")

    # valA: bare magnet at s=0, no coils/shell -> analytic check
    for tag, shell, coils, abc, msc in [("valA", False, False, 60.0, 1.0),
                                        ("valB", True, True, 120.0, 1.0),
                                        ("valC", True, True, 60.0, 0.5)]:
        segs = set()
        emit_doc(w, segs, 0.0, with_shell=shell, with_coils=coils, abc=abc,
                 mesh_scale=msc)
        w(f'mark("{tag} built")')
        w(f'mi_saveas("{z_path(BASE + ".fem")}")')
        if coils:
            for c in M.coils:
                w(f'mi_modifycircprop("{c.name}", 1, 0)')
        w(f'mark("{tag} analyze")')
        w("mi_analyze()")
        w('mark("solved")')
        w("mi_loadsolution()")
        emit_measure(w, tag, 0.0, 0.0)
        w("mo_close()")
        w("mi_close()")

    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def analytic_bare_bz(zs: float) -> float:
    """On-axis Bz of the bare magnet at s=0 (solenoid-sheet, K = Hc)."""
    za, zb, r = -MAG.length / 2, MAG.length / 2, MAG.radius
    c = lambda u: u / math.sqrt(u * u + r * r)
    return (MAG.br / (2 * MAG.mu_r)) * (c(zs - za) - c(zs - zb))


def run():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    csvf = BASE + ".csv"
    for f in (csvf, csvf + ".log"):
        if os.path.exists(f):
            os.remove(f)
    with open(BASE + ".lua", "w") as f:
        f.write(gen_lua())
    subprocess.run(["femm-lua", BASE + ".lua"], check=True, timeout=1500,
                   cwd=PROJECT_ROOT, capture_output=True)
    analyze()


def analyze():
    with open(BASE + ".csv") as f:
        rows = [{k: (v if k == "tag" else float(v))
                 for k, v in r.items()} for r in csv.DictReader(f)]

    def get(tag, s, d):
        return {r["zs"]: r for r in rows
                if r["tag"] == tag and r["s"] == s and r["drive"] == d}

    mt = lambda t: t * 1e3  # T -> mT

    print("== valA: bare magnet vs analytic on-axis formula ==")
    for zs, r in sorted(get("valA", 0, 0).items()):
        an = analytic_bare_bz(zs)
        print(f"  zs={zs:+6.1f}  femm={mt(r['avg_axis']):8.3f} mT  "
              f"analytic={mt(an):8.3f} mT  err={100*(r['avg_axis']/an-1):+5.2f}%")

    ref = get("main", 0.0, 0.0)
    for tag, nm in [("valB", "ABC 60->120"), ("valC", "mesh x0.5")]:
        v = get(tag, 0, 0)
        errs = [100 * (v[zs]["avg_axis"] / ref[zs]["avg_axis"] - 1)
                for zs in SENSOR_Z if abs(ref[zs]["avg_axis"]) > 1e-6]
        nerrs = [100 * (v[zs]["avg_n1"] / ref[zs]["avg_n1"] - 1)
                 for zs in SENSOR_Z if abs(ref[zs]["avg_n1"]) > 1e-7]
        print(f"== {nm}: axis err max {max(abs(e) for e in errs):.2f}%  "
              f"n1 err max {max(abs(e) for e in nerrs):.2f}%")

    print()
    print("== st_w10 sensor candidates (N52 sim scale; xN45SH "
          f"{N45SH_SCALE:.3f} for the real magnet) ==")
    print(f"{'zs':>6} | {'Bz range over stroke (mT)':>26} | "
          f"{'grad mT/mm':>18} | {'lin dev':>7} | {'coil@s=1':>9} | "
          f"{'coil err':>8} | {'nbr dB':>7} | {'nbr err':>8}")
    for zs in SENSOR_Z:
        b = [get("main", s, 0.0)[zs]["avg_axis"] for s in MOVERS]
        # gradient via finite differences over the mover grid
        g = [(b[i + 1] - b[i]) / (MOVERS[i + 1] - MOVERS[i])
             for i in range(len(b) - 1)]
        gmin = min(abs(x) for x in g)
        gmax = max(abs(x) for x in g)
        mono = all(x < 0 for x in g) or all(x > 0 for x in g)
        # linearity: max deviation from the straight line through the ends
        lin = [b[0] + (b[-1] - b[0]) * (s - MOVERS[0]) /
               (MOVERS[-1] - MOVERS[0]) for s in MOVERS]
        span = max(b) - min(b)
        ldev = max(abs(x - y) for x, y in zip(b, lin)) / span if span else 0
        # coil contamination per unit drive (mean of +1/-1 magnitudes, worst s)
        cw = max(0.5 * (abs(get("main", s, 1.0)[zs]["avg_axis"] -
                            get("main", s, 0.0)[zs]["avg_axis"]) +
                        abs(get("main", s, -1.0)[zs]["avg_axis"] -
                            get("main", s, 0.0)[zs]["avg_axis"]))
                 for s in MOVERS)
        # neighbor: stroke-varying part of the unshielded victim field
        nb = [get("main", s, 0.0)[zs]["avg_n1"] for s in MOVERS]
        ndb = max(nb) - min(nb)
        gm = gmin if gmin > 0 else float("nan")
        print(f"{zs:+6.1f} | {mt(min(b)):+9.2f} .. {mt(max(b)):+9.2f}"
              f"{'' if mono else ' NONMONO':>8} | "
              f"{mt(gmin):7.3f} .. {mt(gmax):7.3f} | {100*ldev:6.1f}% | "
              f"{mt(cw):6.3f} mT | {1e3*mt(cw)/mt(gm):5.0f} um | "
              f"{1e3*mt(ndb):5.1f} uT | {1e3*mt(ndb)/mt(gm):5.1f} um")
    print()
    print("coil err = worst apparent position shift at drive s=1 (300 A-t/coil),")
    print("  known-current -> feedforwardable; scale x~3.2 for bus-stall worst case.")
    print("nbr err = apparent shift from ONE nearest neighbor's FULL stroke,")
    print("  unshielded victim upper bound at 19.05 mm pitch.")

    print()
    print("== crosstalk detail at the victim sensor (unshielded, per ring) ==")
    print("ring radii: n1=19.05 (x4: 2 in-row, 2 cross-row), n2=26.94 (x4 diag),")
    print("  n3=38.10 (x4). Adversarial = sum of |dB| over all 12; coil sum =")
    print("  12 neighbors at drive 1 simultaneously (firmware knows currents).")
    print(f"{'zs':>6} | {'dB_mag n1/n2/n3 (uT)':>28} | {'adv err':>8} | "
          f"{'dB_coil n1/n2/n3 (uT)':>28} | {'coil12 err':>10}")
    for zs in SENSOR_Z:
        b = [get("main", s, 0.0)[zs]["avg_axis"] for s in MOVERS]
        g = [(b[i + 1] - b[i]) / (MOVERS[i + 1] - MOVERS[i])
             for i in range(len(b) - 1)]
        gmin = min(abs(x) for x in g)
        dmag, dcoil = [], []
        for nm in ("n1", "n2", "n3"):
            vals = [get("main", s, 0.0)[zs][f"avg_{nm}"] for s in MOVERS]
            dmag.append(max(vals) - min(vals))
            cd = max(0.5 * (abs(get("main", s, 1.0)[zs][f"avg_{nm}"] -
                                get("main", s, 0.0)[zs][f"avg_{nm}"]) +
                            abs(get("main", s, -1.0)[zs][f"avg_{nm}"] -
                                get("main", s, 0.0)[zs][f"avg_{nm}"]))
                     for s in MOVERS)
            dcoil.append(cd)
        adv = 4 * sum(dmag)
        c12 = 4 * sum(dcoil)
        print(f"{zs:+6.1f} | " +
              "/".join(f"{1e6*v:7.1f}" for v in dmag) + " | " +
              f"{1e3*adv/gmin:6.1f} um | " +
              "/".join(f"{1e6*v:7.1f}" for v in dcoil) + " | " +
              f"{1e3*c12/gmin:7.1f} um")


if __name__ == "__main__":
    {"run": run, "analyze": analyze}[sys.argv[1] if len(sys.argv) > 1 else "run"]()
