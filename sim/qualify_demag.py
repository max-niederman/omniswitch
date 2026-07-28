"""Demagnetization qualification + magnet-grade selection for the three
finalist actuators (opp24c3, w_o3_m12_r, st_w03).

Method
------
Each magnet is subdivided into 3 axial slabs (two 1.5 mm end slabs + body)
x 2 radial zones (core r<2.8, rim 2.8<r<4) as separate FEMM blocks sharing
material + magdir; every sub-block gets a unique group id so its volume-
averaged Bz = mo_blockintegral(9)/mo_blockintegral(10) can be read per
block (point evaluation hangs under Wine; block integrals are the
sanctioned path). Sub-block rectangles are emitted on a shared grid with
python-side segment dedup (same trick as MoverModel), so butted faces
(magnet<->washer in w_o3_m12_r) share identical split segments and no
coincident-collinear segments arise -> no boundary inset needed.

For temp T in {60,80,100,120} C and grade G in {N52,N45SH,N38UH} the magnet
material Br is set to GRADES[G].br(T) (Hc = Br/(mu_r*mu0)), and the model
is solved at mover z in {0, +2.5} (z=-2.5 follows from the z -> -z +
magnet-swap + s -> -s symmetry of all three designs) with coil scale
s in {-s_peak, 0, +s_peak}, where s_peak = s_hold_worst / 0.30 = the
bus-limited stall overdrive implied by analyze.py's 30% hold duty.
Solving both signs of s covers "the coil polarity that reverses field in
the magnet" for BOTH magnets of the opposed pairs; s=0 gives the
permanent (unpowered) mutual-reverse-field state.

Margin: H_d = (Bz_local_along_M - Br(T)) / (mu0 mu_r); a 1.5x local
peaking factor is applied on |H_d| for corner/edge hotspots that a
6-block average cannot resolve: margin_peak = Hk(T) / (1.5 |H_d|).
margin_avg (no peaking) is reported alongside.

Fz on the whole mover (all magnet sub-block groups + washer group) is also
recorded per solve -> empirical force gain per (grade, temp) -> hold power
of the chosen grade incl. steel nonlinearity, cross-checked against the
analytic power_scale() ~ 1/Br(T)^2.

Retention: static WST z-forces at s=0, 20 C, N52 (Br 1.445 T per GRADES;
the finalist sweeps used 1.43 T -> forces here read ~2% high vs those):
per-magnet forces for both opposed pairs (each magnet its own group), and
for w_o3_m12_r the magnet->washer decomposition. w_o3 retention models get
an explicit 0.1 mm adhesive bond-line air gap between magnet and washer
(real glue line; also gives the WST group an air wrap -> valid force
partitioning). Gap regions are enclosed by a small segment lip and meshed
at 0.05 mm.

s_peak derivations (from the frozen finalist fine runs):
  opp24c3    : worst gain 5.970 N/unit (opp24c3f_acd4, z=+/-2.5)
               -> s_hold 0.1340, s_peak 0.4467  (mid 268 A-t, outers 134)
  w_o3_m12_r : worst gain 2.8409 N/unit (w_o3_m12_r_wcgF, z=+2.5)
               -> s_hold 0.2816, s_peak 0.9387  (mid 282 A-t, outers 108)
  st_w03     : worst s(+0.8 N) 0.3790 (st_w03_fine_stloptfin, z=+2.5)
               -> s_peak 1.2633  (379 A-t per coil)

Run (from repo root, inside nix develop, ALWAYS under an outer timeout):
    python sim/qualify_demag.py smoke <design>      # 1 condition sanity
    python sim/qualify_demag.py run <design>        # full T x G matrix
    python sim/qualify_demag.py analyze             # tables + margins CSV
    python sim/qualify_demag.py retention           # all retention models
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from dataclasses import dataclass, field

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import MU0, PROJECT_ROOT, RESULTS_DIR, Coil, z_path  # noqa: E402
from magnets import GRADES, power_scale  # noqa: E402

TEMPS = [60.0, 80.0, 100.0, 120.0]
GRADE_LIST = ["N52", "N45SH", "N38UH"]
Z_LIST = [0.0, 2.5]
CORE_R = 2.8        # core/rim split radius
END_SLAB = 1.5      # end-slab thickness (mm)
MESH_MAG = 0.25     # mesh size in magnet sub-blocks (mm)
PEAKING = 1.5       # local |H_d| peaking factor for corner hotspots
F_TARGET = 0.8
HOLD_DUTY = 0.30
MARGIN_REQ = 1.3    # required margin at 100 C


@dataclass
class Part:
    name: str          # "A"/"B" (A = lower magnet)
    z_off: float       # part center relative to mover reference z
    length: float
    dirsign: int       # +1 = magnetized +z (magdir 90), -1 = -z (magdir 270)
    facing: int        # +1: facing (opposed) end is top; -1: bottom; 0: none
    radius: float = 4.0


@dataclass
class Design:
    name: str
    parts: list
    coils: list
    washer: tuple | None = None   # (z_bot, z_top) mover steel r 0..4
    shell: tuple | None = None    # (r_in, r_out, z_bot, z_top) stator steel
    s_peak: float = 0.0
    s_hold: float = 0.0
    abc_radius: float = 60.0

    def mover_groups(self):
        gs = []
        for mi, _ in enumerate(self.parts):
            gs += [10 * (mi + 1) + k for k in range(6)]
        if self.washer:
            gs.append(30)
        return gs


DESIGNS = {
    "opp24c3": Design(
        name="opp24c3",
        parts=[Part("A", -6.1, 12.0, +1, facing=+1),
               Part("B", +6.1, 12.0, -1, facing=-1)],
        coils=[Coil("lo", 4.5, 7.75, -16.5, -8.0, ni=-300.0),
               Coil("mid", 4.5, 7.75, -6.5, 6.5, ni=600.0),
               Coil("hi", 4.5, 7.75, 8.0, 16.5, ni=-300.0)],
        s_hold=0.1340, s_peak=0.4467),
    "w_o3_m12_r": Design(
        name="w_o3_m12_r",
        parts=[Part("A", -7.0, 12.0, +1, facing=+1),
               Part("B", +7.0, 12.0, -1, facing=-1)],
        coils=[Coil("mid", 4.75, 7.75, -5.0, 5.0, ni=300.0),
               Coil("hi", 4.75, 7.75, 5.5, 15.5, ni=-115.0),
               Coil("lo", 4.75, 7.75, -15.5, -5.5, ni=-115.0)],
        washer=(-1.0, 1.0),
        s_hold=0.2816, s_peak=0.9387),
    "st_w03": Design(
        name="st_w03",
        parts=[Part("A", 0.0, 12.0, +1, facing=0)],
        coils=[Coil("lo", 4.75, 7.65, -12.5, -0.5, ni=-300.0),
               Coil("hi", 4.75, 7.65, 0.5, 12.5, ni=+300.0)],
        shell=(7.7, 8.0, -13.0, 13.0),
        s_hold=0.3790, s_peak=1.2633),
}


def slab_name(part: Part, slab_i: int) -> str:
    # slab_i: 0 = bottom 1.5 mm, 1 = body, 2 = top 1.5 mm
    if part.facing == 0:
        return ["bot", "body", "top"][slab_i]
    if part.facing > 0:
        return ["far", "body", "facing"][slab_i]
    return ["facing", "body", "far"][slab_i]


def subblocks(part: Part):
    """Yield (slab_i, zone_i, r1, r2, z1, z2) in mover-relative coords."""
    zb = part.z_off - part.length / 2
    zt = part.z_off + part.length / 2
    slabs = [(zb, zb + END_SLAB), (zb + END_SLAB, zt - END_SLAB),
             (zt - END_SLAB, zt)]
    zones = [(0.0, CORE_R), (CORE_R, part.radius)]
    for si, (z1, z2) in enumerate(slabs):
        for zi, (r1, r2) in enumerate(zones):
            yield si, zi, r1, r2, z1, z2


# ------------------------------------------------------------------ lua emit

def _mk_writers(L: list):
    w = L.append
    segs: set = set()

    def rect(r1, z1, r2, z2):
        for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                             (r2, z2, r1, z2), (r1, z2, r1, z1)]:
            key = tuple(sorted([(round(a, 6), round(b, 6)),
                                (round(c, 6), round(d, 6))]))
            if key in segs:
                continue
            segs.add(key)
            w(f"mi_addnode({a:.6g}, {b:.6g})")
            w(f"mi_addnode({c:.6g}, {d:.6g})")
            w(f"mi_addsegment({a:.6g}, {b:.6g}, {c:.6g}, {d:.6g})")

    def seg(a, b, c, d):
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

    return w, rect, seg, label


def emit_demag_doc(L: list, d: Design, z: float, br: float, mu_r: float):
    w, rect, _seg, label = _mk_writers(L)
    hc = br / (mu_r * MU0)
    w("newdocument(0)")
    w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w('mi_getmaterial("Copper")')
    w(f'mi_addmaterial("Mag", {mu_r}, {mu_r}, {hc:.1f}, '
      "0, 0.667, 0, 0, 1, 0, 0, 0)")
    if d.washer or d.shell:
        w('mi_getmaterial("1018 Steel")')

    for mi, p in enumerate(d.parts):
        magdir = 90 if p.dirsign > 0 else 270
        for si, zi, r1, r2, z1, z2 in subblocks(p):
            gid = 10 * (mi + 1) + si * 2 + zi
            rect(r1, z + z1, r2, z + z2)
            label((r1 + r2) / 2, z + (z1 + z2) / 2, "Mag",
                  magdir=magdir, group=gid, mesh=MESH_MAG)

    if d.washer:
        wz1, wz2 = d.washer
        # split at CORE_R so butted faces share identical grid segments
        for (r1, r2) in [(0.0, CORE_R), (CORE_R, 4.0)]:
            rect(r1, z + wz1, r2, z + wz2)
            label((r1 + r2) / 2, z + (wz1 + wz2) / 2, "1018 Steel",
                  group=30, mesh=0.3)

    for c in d.coils:
        rect(c.r_in, c.z_bot, c.r_out, c.z_top)
        w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
        label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
              circuit=c.name, turns=1)

    if d.shell:
        sr1, sr2, sz1, sz2 = d.shell
        rect(sr1, sz1, sr2, sz2)
        label((sr1 + sr2) / 2, (sz1 + sz2) / 2, "1018 Steel", group=0,
              mesh=0.3)

    w(f"mi_makeABC(7, {d.abc_radius:.6g}, 0, 0, 0)")
    label(d.abc_radius * 0.5, d.abc_radius * 0.5, "Air")


def demag_lua(d: Design, grades, temps, zs, blocks_csv, force_csv,
              fem_file) -> str:
    L: list[str] = []
    w = L.append
    w(f'LOG = "{z_path(blocks_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'hb = openfile("{z_path(blocks_csv)}", "w")')
    w('write(hb, "z,grade,temp,s,part,slab,zone,dir,group,intBz,vol\\n")')
    w(f'hf = openfile("{z_path(force_csv)}", "w")')
    w('write(hf, "z,grade,temp,s,Fz\\n")')

    for z in zs:
        for gname in grades:
            g = GRADES[gname]
            for T in temps:
                br = g.br(T)
                emit_demag_doc(L, d, z, br, g.mu_r)
                w(f'mark("z={z:.6g} {gname} T={T:.6g} built")')
                w(f'mi_saveas("{z_path(fem_file)}")')
                for s in (-d.s_peak, 0.0, d.s_peak):
                    for c in d.coils:
                        w(f'mi_modifycircprop("{c.name}", 1, '
                          f"{c.ni * s:.6g})")
                    w(f'mark("z={z:.6g} {gname} T={T:.6g} s={s:.6g} '
                      'analyze")')
                    w("mi_analyze()")
                    w('mark("solved")')
                    w("mi_loadsolution()")
                    for grp in d.mover_groups():
                        w(f"mo_groupselectblock({grp})")
                    w("fz = mo_blockintegral(19)")
                    w("mo_clearblock()")
                    w(f'write(hf, "{z:.6g},{gname},{T:.6g},{s:.6g},", '
                      'fz, "\\n")')
                    for mi, p in enumerate(d.parts):
                        for si, zi, _r1, _r2, _z1, _z2 in subblocks(p):
                            gid = 10 * (mi + 1) + si * 2 + zi
                            zone = "core" if zi == 0 else "rim"
                            w(f"mo_groupselectblock({gid})")
                            w("ib = mo_blockintegral(9)")
                            w("vol = mo_blockintegral(10)")
                            w("mo_clearblock()")
                            w(f'write(hb, "{z:.6g},{gname},{T:.6g},'
                              f"{s:.6g},{p.name},{slab_name(p, si)},"
                              f'{zone},{p.dirsign},{gid},", ib, ",", '
                              'vol, "\\n")')
                    w("mo_close()")
                w("mi_close()")
    w("closefile(hb)")
    w("closefile(hf)")
    w("quit()")
    return "\n".join(L) + "\n"


def run_lua(lua_text: str, base: str, timeout: int = 3000):
    lua_file = base + ".lua"
    with open(lua_file, "w") as f:
        f.write(lua_text)
    subprocess.run(["femm-lua", lua_file], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)


def run_design(d: Design, grades=None, temps=None, zs=None, tag="demag"):
    grades = grades or GRADE_LIST
    temps = temps or TEMPS
    zs = zs or Z_LIST
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, f"{d.name}_{tag}")
    blocks_csv, force_csv = base + "_blocks.csv", base + "_force.csv"
    for p in (blocks_csv, force_csv):
        if os.path.exists(p):
            os.remove(p)
    run_lua(demag_lua(d, grades, temps, zs, blocks_csv, force_csv,
                      base + ".fem"), base)
    return blocks_csv, force_csv


# ------------------------------------------------------------------ analysis

def margin_rows(blocks_csv: str, design: str):
    out = []
    with open(blocks_csv) as f:
        for r in csv.DictReader(f):
            g = GRADES[r["grade"]]
            T = float(r["temp"])
            avg_bz = float(r["intBz"]) / float(r["vol"])
            bz_local = float(r["dir"]) * avg_bz   # along magnetization
            h_d = (bz_local - g.br(T)) / (MU0 * g.mu_r)
            if h_d >= 0:
                m_avg = m_pk = float("inf")
            else:
                m_avg = g.hk(T) / -h_d
                m_pk = g.hk(T) / (PEAKING * -h_d)
            out.append(dict(design=design, grade=r["grade"], temp=T,
                            z=float(r["z"]), s=float(r["s"]),
                            part=r["part"], slab=r["slab"], zone=r["zone"],
                            avg_bz=bz_local, margin_avg=m_avg,
                            margin_peak=m_pk))
    return out


def analyze():
    all_rows = []
    force = {}
    for name in DESIGNS:
        base = os.path.join(RESULTS_DIR, f"{name}_demag")
        all_rows += margin_rows(base + "_blocks.csv", name)
        with open(base + "_force.csv") as f:
            for r in csv.DictReader(f):
                force[(name, r["grade"], float(r["temp"]), float(r["z"]),
                       float(r["s"]))] = float(r["Fz"])

    # min margin per (design, grade, temp)
    mins = {}
    for r in all_rows:
        k = (r["design"], r["grade"], r["temp"])
        if k not in mins or r["margin_peak"] < mins[k]["margin_peak"]:
            mins[k] = r
    # min margin at s=0 only (unpowered mutual-field state)
    mins0 = {}
    for r in all_rows:
        if r["s"] != 0.0:
            continue
        k = (r["design"], r["grade"], r["temp"])
        if k not in mins0 or r["margin_peak"] < mins0[k]["margin_peak"]:
            mins0[k] = r

    out_csv = os.path.join(RESULTS_DIR, "demag_margins.csv")
    with open(out_csv, "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["design", "grade", "temp", "margin_peak_min",
                       "margin_avg_min", "margin_peak_s0", "worst_z",
                       "worst_s", "worst_part", "worst_slab", "worst_zone",
                       "worst_avgBz"])
        for name in DESIGNS:
            for gname in GRADE_LIST:
                for T in TEMPS:
                    r = mins[(name, gname, T)]
                    r0 = mins0[(name, gname, T)]
                    wcsv.writerow([name, gname, T,
                                   f"{r['margin_peak']:.3f}",
                                   f"{r['margin_avg']:.3f}",
                                   f"{r0['margin_peak']:.3f}",
                                   r["z"], r["s"], r["part"], r["slab"],
                                   r["zone"], f"{r['avg_bz']:.4f}"])

    print("=== min demag margin (1.5x peaking on |H_d|) per design/grade/"
          "temp ===")
    print("margin_avg in parens; worst sub-block noted; s0 = unpowered")
    for name in DESIGNS:
        print(f"\n-- {name} (s_peak={DESIGNS[name].s_peak}) --")
        hdr = "grade    " + "".join(f"  T={T:<5.0f}       " for T in TEMPS)
        print(hdr)
        for gname in GRADE_LIST:
            cells = []
            for T in TEMPS:
                r = mins[(name, gname, T)]
                cells.append(f"  {r['margin_peak']:5.2f}({r['margin_avg']:5.2f})")
            r100 = mins[(name, gname, 100.0)]
            print(f"{gname:8s}" + "".join(cells)
                  + f"   worst@100C: {r100['part']}/{r100['slab']}/"
                    f"{r100['zone']} z={r100['z']} s={r100['s']:+.3f} "
                    f"avgBz={r100['avg_bz']:+.3f} T")
        for gname in GRADE_LIST:
            r0 = mins0[(name, gname, 100.0)]
            print(f"   {gname:8s} unpowered(s=0)@100C: "
                  f"margin_peak={r0['margin_peak']:5.2f} "
                  f"({r0['part']}/{r0['slab']}/{r0['zone']}, "
                  f"avgBz={r0['avg_bz']:+.3f} T)")

    print("\n=== (design, grade) with margin_peak >= "
          f"{MARGIN_REQ} at 100 C ===")
    for name in DESIGNS:
        ok = [g for g in GRADE_LIST
              if mins[(name, g, 100.0)]["margin_peak"] >= MARGIN_REQ]
        print(f"{name:12s}: {ok if ok else 'NONE'}")

    # empirical hold power per grade/temp from the recorded Fz gains
    print("\n=== hold power (0.8 N worst dir, from FEMM gains at z in "
          "{0, 2.5}) ===")
    print("P[W] at fill 0.6 (fill 0.5 = x1.2); analytic power_scale vs "
          "N52@same-T in parens")
    for name, d in DESIGNS.items():
        print(f"\n-- {name} --")
        for gname in GRADE_LIST:
            cells = []
            for T in TEMPS:
                gains = []
                for z in Z_LIST:
                    fp = force[(name, gname, T, z, d.s_peak)]
                    fm = force[(name, gname, T, z, -d.s_peak)]
                    gains.append((fp - fm) / (2 * d.s_peak))
                gmin = min(gains)
                s_hold = F_TARGET / gmin
                p06 = sum(c.power(c.ni * s_hold) for c in d.coils)
                ps = power_scale("N52", gname, T)
                cells.append(f"  {p06:5.3f}({ps:4.2f})")
            print(f"{gname:8s}" + "".join(
                f" T={T:<3.0f}:{c}" for T, c in zip(TEMPS, cells)))
    print(f"\nsummary CSV: {out_csv}")


# ------------------------------------------------------------------ retention

RET_BR = GRADES["N52"].br20   # 1.445 T (finalist sweeps used 1.43 -> ~2%)
RET_MU = GRADES["N52"].mu_r


def emit_ret_doc(L, parts, washer, gaps, abc=60.0):
    """parts: (name, z_off, length, dirsign, group); washer: (z1, z2, group)
    or None; gaps: list of (z1, z2) air gaps to enclose+fine-mesh."""
    w, rect, seg, label = _mk_writers(L)
    hc = RET_BR / (RET_MU * MU0)
    w("newdocument(0)")
    w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w(f'mi_addmaterial("Mag", {RET_MU}, {RET_MU}, {hc:.1f}, '
      "0, 0.667, 0, 0, 1, 0, 0, 0)")
    if washer:
        w('mi_getmaterial("1018 Steel")')
    for (name, z_off, ln, dirsign, grp) in parts:
        rect(0.0, z_off - ln / 2, 4.0, z_off + ln / 2)
        label(2.0, z_off, "Mag", magdir=90 if dirsign > 0 else 270,
              group=grp, mesh=0.2)
    if washer:
        wz1, wz2, grp = washer
        rect(0.0, wz1, 4.0, wz2)
        label(2.0, (wz1 + wz2) / 2, "1018 Steel", group=grp, mesh=0.15)
    for (gz1, gz2) in gaps:
        # lip enclosing the gap region so it can carry a 0.05 mm mesh label
        seg(4.0, gz1, 4.3, gz1)
        seg(4.3, gz1, 4.3, gz2)
        seg(4.3, gz2, 4.0, gz2)
        label(2.0, (gz1 + gz2) / 2, "Air", mesh=0.05)
    w(f"mi_makeABC(7, {abc:.6g}, 0, 0, 0)")
    label(abc * 0.5, abc * 0.5, "Air")


RET_MODELS = {
    # opp24c3 mover at rest: two D8x12 N-to-N across the 0.2 mm spacer
    "opp24c3_pair": dict(
        parts=[("A", -6.1, 12.0, +1, 1), ("B", +6.1, 12.0, -1, 2)],
        washer=None, gaps=[(-0.1, 0.1)],
        groups={1: "magnet_A(lower)", 2: "magnet_B(upper)"}),
    # w_o3_m12_r full mover, 0.1 mm bond-line gaps each side of the washer
    "w_o3_full": dict(
        parts=[("A", -7.1, 12.0, +1, 1), ("B", +7.1, 12.0, -1, 2)],
        washer=(-1.0, 1.0, 3), gaps=[(-1.1, -1.0), (1.0, 1.1)],
        groups={1: "magnet_A(lower)", 2: "magnet_B(upper)", 3: "washer"}),
    # assembly step 1: one magnet approaching the bare washer
    "w_o3_soloA": dict(
        parts=[("A", -7.1, 12.0, +1, 1)],
        washer=(-1.0, 1.0, 3), gaps=[(-1.1, -1.0)],
        groups={1: "magnet_A(lower)", 3: "washer"}),
    # decomposition: the two magnets at assembled spacing, washer removed
    "w_o3_magsonly": dict(
        parts=[("A", -7.1, 12.0, +1, 1), ("B", +7.1, 12.0, -1, 2)],
        washer=None, gaps=[],
        groups={1: "magnet_A(lower)", 2: "magnet_B(upper)"}),
}


def retention():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, "retention_qd")
    out_csv = base + ".csv"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    L: list[str] = []
    w = L.append
    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'h = openfile("{z_path(out_csv)}", "w")')
    w('write(h, "model,group,Fz\\n")')
    for name, spec in RET_MODELS.items():
        emit_ret_doc(L, spec["parts"], spec["washer"], spec["gaps"])
        w(f'mark("{name} built")')
        w(f'mi_saveas("{z_path(base)}_{name}.fem")')
        w(f'mark("{name} analyze")')
        w("mi_analyze()")
        w('mark("solved")')
        w("mi_loadsolution()")
        for grp in spec["groups"]:
            w(f"mo_groupselectblock({grp})")
            w("fz = mo_blockintegral(19)")
            w("mo_clearblock()")
            w(f'write(h, "{name},{grp},", fz, "\\n")')
        w("mo_close()")
        w("mi_close()")
    w("closefile(h)")
    w("quit()")
    run_lua("\n".join(L) + "\n", base, timeout=1800)

    print(f"=== retention forces (s=0, 20 C, N52 Br={RET_BR} T) ===")
    print("(+Fz = up; A is the lower magnet; w_o3 models have a 0.1 mm "
          "bond-line gap)")
    with open(out_csv) as f:
        for r in csv.DictReader(f):
            gname = RET_MODELS[r["model"]]["groups"][int(r["group"])]
            print(f"{r['model']:14s} {gname:16s} Fz = "
                  f"{float(r['Fz']):+8.2f} N")
    print(f"CSV: {out_csv}")


# ------------------------------------------------------------------ cli

def main(argv):
    mode = argv[0]
    if mode == "smoke":
        d = DESIGNS[argv[1]]
        b, fcsv = run_design(d, grades=["N52"], temps=[80.0], zs=[0.0],
                             tag="demagsmoke")
        for row in margin_rows(b, d.name):
            print(row)
        with open(fcsv) as f:
            print(f.read())
    elif mode == "run":
        d = DESIGNS[argv[1]]
        run_design(d)
        print(f"done {d.name}")
    elif mode == "analyze":
        analyze()
    elif mode == "retention":
        retention()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1:])
