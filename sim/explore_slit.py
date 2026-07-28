"""Axial-slit trade study for the frozen st_w10 shell (see docs/DESIGN.md).

Question: full-length axial slit(s) in the 1.0 mm 1018 shell to kill eddy
currents — single slit, two opposed slits, and end-bridged slits.

Part AC (harmonic, extends sim/explore_steel_ac.py to PWM frequencies):
  st_w10 geometry, LINEAR steel (mur 325 = measured differential mu at the
  1.21 T wall bias; mur 529 control lives in results/st_w10_stlac.csv),
  magnet linear mur 1.05 Hc=0 sigma 0.667 MS/m, coils at base NI=300.
  freqs 0..50 kHz, steel mesh sized to the skin depth (delta/3, floor 0.03).
  Configs:
    ctrl_air  : no steel, sigma_mag=0 (method control)
    mag_eddy  : no steel, magnet eddies only  = the no-shell R(f) baseline
    steel325  : full shell, unslit
    slit325   : + parallel 0 A circuit on the shell = net circumferential
                current forced to zero (the topological effect of a slit on
                the NET winding mode; z-local zero-net-winding eddy loops
                remain — see analysis: for the series-opposed coil drive the
                induced E_phi is odd in z and the net mode is NOT excited,
                so slit325 == steel325 is the expected/verified outcome)
    bridge325 : only the 1.5 mm end rings conductive (sigma 0 elsewhere) =
                axisymmetric bound of the residual shorted-turn modes a
                bridged (slit-stopped-short-of-ends) tube keeps as REAL
                complete rings
    steel325cm/slit325cm : both coils driven SAME sign (common mode) — the
                excitation a slit actually kills; demonstrates the contrast.
  Rows are appended one solve at a time to results/st_w10_slitac.csv and
  existing (config,freq) rows are skipped -> resumable after a Wine hang.

Part SIDE (DC, magnetostatic): B_n(z) on the shell inner surface via 26
  contour windows of 1 mm at r=6.975 (mid of the 0.05 mm coil-shell air gap),
  mo_lineintegral(0) for flux and (5) for the exact (B.n)^2 window average
  (point evaluation hangs under Wine). Mover z in {-2.5, 0, +2.5}, drive
  s in {0, +/-0.3809} (worst-direction hold from converge_shell report).
  Side-load of a slit of kerf w = Maxwell tension of the missing wedge:
  F = sum avg(Bn^2)_i * dz * w / (2 mu0)  (worst case: full stress removal).

Run inside `nix develop` from the repo root, ALWAYS under `timeout`:
    python sim/explore_slit.py side            # 9 solves, ~5 min
    python sim/explore_slit.py side-analyze
    python sim/explore_slit.py ac-all          # ~40 solves, resumable
    python sim/explore_slit.py ac-analyze
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from femm import MU0, PROJECT_ROOT, RESULTS_DIR, z_path  # noqa: E402
from explore_steel_shell import VARIANTS  # noqa: E402

M = VARIANTS["st_w10"]
S_HOLD = 0.3809          # worst-direction 0.8 N hold scale (converge report)
MUR_BIAS = 325.0         # differential mu_r of 1018 at the 1.21 T wall bias
SIG_STEEL = 5.8          # MS/m
SIG_MAG = 0.667          # MS/m sintered NdFeB
BRIDGE = 1.5             # mm, axial length of each end bridge (2c)
SHELL = M.steels[0]      # r 7.0-8.0, z +/-13

AC_CSV = os.path.join(RESULTS_DIR, "st_w10_slitac.csv")
SIDE_CSV = os.path.join(RESULTS_DIR, "st_w10_slitside.csv")

FREQS = [0.0, 20.0, 100.0, 500.0, 1000.0, 5000.0, 25000.0, 50000.0]
FREQS_CM = [0.0, 100.0, 500.0, 25000.0]

# config -> (with_steel, mag_sigma, slit_circuit, bridge_only, common_mode)
AC_CONFIGS = {
    "ctrl_air":   (False, 0.0,      False, False, False),
    "mag_eddy":   (False, SIG_MAG,  False, False, False),
    "steel325":   (True,  SIG_MAG,  False, False, False),
    "slit325":    (True,  SIG_MAG,  True,  False, False),
    "bridge325":  (True,  SIG_MAG,  False, True,  False),
    "steel325cm": (True,  SIG_MAG,  False, False, True),
    "slit325cm":  (True,  SIG_MAG,  True,  False, True),
}
AC_PLAN = [(c, f) for c in ("ctrl_air", "mag_eddy", "steel325", "slit325",
                            "bridge325") for f in FREQS] + \
          [(c, f) for c in ("steel325cm", "slit325cm") for f in FREQS_CM]

AC_COLS = ("config,freq,absI_lo,argI_lo,absV_lo,argV_lo,absLam_lo,argLam_lo,"
           "absI_hi,argI_hi,absV_hi,argV_hi,absLam_hi,argLam_hi,"
           "Psteel,absIsteel,Pmag")


def skin_mm(f: float, mur: float = MUR_BIAS) -> float:
    if f <= 0:
        return 1e9
    return 1e3 / math.sqrt(math.pi * f * mur * 4e-7 * math.pi
                           * SIG_STEEL * 1e6)


def steel_mesh(f: float) -> float:
    return min(0.2, max(0.03, skin_mm(f) / 3.0))


class Geo:
    """Segment-deduping emitter (bridge slabs butt against the middle)."""

    def __init__(self, w):
        self.w = w
        self.segs = set()

    def seg(self, r1, z1, r2, z2):
        key = tuple(sorted([(round(r1, 6), round(z1, 6)),
                            (round(r2, 6), round(z2, 6))]))
        if key in self.segs:
            return
        self.segs.add(key)
        self.w(f"mi_addnode({r1:.6g}, {z1:.6g})")
        self.w(f"mi_addnode({r2:.6g}, {z2:.6g})")
        self.w(f"mi_addsegment({r1:.6g}, {z1:.6g}, {r2:.6g}, {z2:.6g})")

    def rect(self, r1, z1, r2, z2):
        for (a, b, c, d) in [(r1, z1, r2, z1), (r2, z1, r2, z2),
                             (r2, z2, r1, z2), (r1, z2, r1, z1)]:
            self.seg(a, b, c, d)

    def label(self, r, z, mat, circuit="", magdir=0, group=0, turns=0,
              mesh=0.0):
        automesh = 1 if mesh <= 0 else 0
        self.w(f"mi_addblocklabel({r:.6g}, {z:.6g})")
        self.w(f"mi_selectlabel({r:.6g}, {z:.6g})")
        self.w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, '
               f'"{circuit}", {magdir}, {group}, {turns})')
        self.w("mi_clearselected()")


# ------------------------------------------------------------------ AC part
def ac_one_lua(cfg: str, f: float, out_csv: str, fem: str) -> str:
    with_steel, msig, slit, bridge, cm = AC_CONFIGS[cfg]
    mag = M.magnet
    L: list[str] = []
    w = L.append
    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')

    w("newdocument(0)")
    w(f'mi_probdef({f:.6g}, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w('mi_addmaterial("CuAC", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)')
    w(f'mi_addmaterial("MagnetAC", {mag.mu_r}, {mag.mu_r}, 0, 0, '
      f"{msig:.6g}, 0, 0, 1, 0, 0, 0)")
    if with_steel:
        w(f'mi_addmaterial("SteelC", {MUR_BIAS:.6g}, {MUR_BIAS:.6g}, 0, 0, '
          f"{SIG_STEEL:.6g}, 0, 0, 1, 0, 0, 0)")
        w(f'mi_addmaterial("SteelNC", {MUR_BIAS:.6g}, {MUR_BIAS:.6g}, 0, 0, '
          "0, 0, 0, 1, 0, 0, 0)")
    steel_circ = ""
    if with_steel and slit:
        w('mi_addcircprop("eddy0", 0, 0)')      # parallel, 0 A net
        steel_circ = "eddy0"

    g = Geo(w)
    g.rect(0, -mag.length / 2, mag.radius, mag.length / 2)
    g.label(mag.radius / 2, 0, "MagnetAC", group=1, mesh=0.4)

    for c in M.coils:
        g.rect(c.r_in, c.z_bot, c.r_out, c.z_top)
        ni = abs(c.ni) if cm else c.ni
        w(f'mi_addcircprop("{c.name}", {ni:.6g}, 1)')
        g.label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "CuAC",
                circuit=c.name, turns=1)

    if with_steel:
        s = SHELL
        msz = steel_mesh(f)
        if bridge:
            zb = s.z_top - BRIDGE
            g.rect(s.r_in, -s.z_top, s.r_out, -zb)
            g.rect(s.r_in, -zb, s.r_out, zb)
            g.rect(s.r_in, zb, s.r_out, s.z_top)
            rm = (s.r_in + s.r_out) / 2
            g.label(rm, -(s.z_top + zb) / 2, "SteelC", circuit=steel_circ,
                    group=2, mesh=msz)
            g.label(rm, (s.z_top + zb) / 2, "SteelC", circuit=steel_circ,
                    group=2, mesh=msz)
            g.label(rm, 0, "SteelNC", group=2, mesh=0.2)
        else:
            g.rect(s.r_in, s.z_bot, s.r_out, s.z_top)
            g.label((s.r_in + s.r_out) / 2, 0, "SteelC",
                    circuit=steel_circ, group=2, mesh=msz)

    w(f"mi_makeABC(7, {M.abc_radius:.6g}, 0, 0, 0)")
    g.label(M.abc_radius * 0.5, M.abc_radius * 0.5, "Air")

    w(f'mark("{cfg} f={f:.6g} built (steel mesh '
      f'{steel_mesh(f):.3g} mm, skin {skin_mm(f):.3g} mm)")')
    w(f'mi_saveas("{z_path(fem)}")')
    w("mi_analyze()")
    w(f'mark("{cfg} f={f:.6g} solved")')
    w("mi_loadsolution()")
    w(f'write(handle, "{cfg}", ",", {f:.6g})')
    for c in M.coils:
        w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
        w('write(handle, ",", abs(ic), ",", arg(ic), ",", abs(vc), ",", '
          'arg(vc), ",", abs(lam), ",", arg(lam))')
    if with_steel:
        w("mo_groupselectblock(2)")
        w("ps = mo_blockintegral(4)")       # resistive (eddy) losses, W
        w("iw = mo_blockintegral(7)")       # net current, A (complex)
        w("mo_clearblock()")
        w('write(handle, ",", abs(ps), ",", abs(iw))')
    else:
        w('write(handle, ",0,0")')
    w("mo_groupselectblock(1)")
    w("pm = mo_blockintegral(4)")
    w("mo_clearblock()")
    w('write(handle, ",", abs(pm), "\\n")')
    w("mo_close()")
    w("mi_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def ac_done() -> set:
    if not os.path.exists(AC_CSV):
        return set()
    with open(AC_CSV) as fh:
        return {(r["config"], float(r["freq"]))
                for r in csv.DictReader(fh)}


def ac_run_one(cfg: str, f: float, timeout: int = 2400):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, "st_w10_slitac_one")
    lua, out1, fem = base + ".lua", base + ".csv", base + ".fem"
    if os.path.exists(out1):
        os.remove(out1)
    with open(lua, "w") as fh:
        fh.write(ac_one_lua(cfg, f, out1, fem))
    subprocess.run(["femm-lua", lua], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    with open(out1) as fh:
        row = fh.read().strip()
    new = not os.path.exists(AC_CSV)
    with open(AC_CSV, "a") as fh:
        if new:
            fh.write(AC_COLS + "\n")
        fh.write(row + "\n")
    print(f"done {cfg} f={f:g}")


def ac_all():
    done = ac_done()
    for (cfg, f) in AC_PLAN:
        if (cfg, f) in done:
            print(f"skip {cfg} f={f:g} (have)")
            continue
        ac_run_one(cfg, f)


def _r_of(row, coil):
    """Re(Z) of one 1-turn circuit, ohms (convention-free)."""
    dv = float(row[f"absV_{coil}"])
    di = float(row[f"absI_{coil}"])
    if di == 0:
        return 0.0
    ph = float(row[f"argV_{coil}"]) - float(row[f"argI_{coil}"])
    return dv / di * math.cos(ph)


def ac_analyze():
    with open(AC_CSV) as fh:
        rows = list(csv.DictReader(fh))
    by = {}
    for r in rows:
        by.setdefault(r["config"], []).append(r)
    for cfg, rr in by.items():
        rr.sort(key=lambda r: float(r["freq"]))
        ref = next(r for r in rr if float(r["freq"]) == 0.0)
        lam0 = float(ref["absLam_hi"])
        arg0 = float(ref["argLam_hi"])
        print(f"-- {cfg} (|lam0|={lam0:.4e} Wb)")
        print("     f[Hz]  |lam|/|lam0|  phase[deg]  tau_eq[us]   "
              "R1t_lo+hi[ohm]  Psteel[W]@300At  |Isteel|[A]  Pmag[W]")
        for r in rr:
            f = float(r["freq"])
            rel = float(r["absLam_hi"]) / lam0
            phi = float(r["argLam_hi"]) - arg0
            while phi > math.pi:
                phi -= 2 * math.pi
            while phi < -math.pi:
                phi += 2 * math.pi
            tau = math.tan(-phi) / (2 * math.pi * f) * 1e6 if f > 0 else 0.0
            rtot = _r_of(r, "lo") + _r_of(r, "hi")
            print(f"  {f:8.0f}  {rel:10.4f}  {math.degrees(phi):+9.3f}  "
                  f"{tau:9.3f}   {rtot:12.4e}   {float(r['Psteel']):11.4e}  "
                  f"{float(r['absIsteel']):9.3f}  {float(r['Pmag']):.3e}")
    # ripple-loss synthesis at 25 kHz
    print("\n== 25 kHz PWM ripple synthesis ==")
    n_turns, i_pp = 480.0, 0.022
    ni_pp = n_turns * i_pp
    a1_pk = 8 / math.pi ** 2 * ni_pp / 2      # triangle fundamental, peak A-t
    a1_rms = a1_pk / math.sqrt(2)
    print(f"22 mA p-p triangular at N={n_turns:.0f} -> NI_pp = {ni_pp:.2f} "
          f"A-t; fundamental peak {a1_pk:.2f} A-t, rms {a1_rms:.2f} A-t "
          "(harmonics: amp 1/n^2 -> <3% extra loss, ignored)")

    def at(cfg, f):
        return next(r for r in by[cfg] if float(r["freq"]) == f)

    for cfg in ("steel325", "slit325", "bridge325"):
        if cfg not in by:
            continue
        try:
            r25 = at(cfg, 25000.0)
            rb = at("mag_eddy", 25000.0)
        except StopIteration:
            continue
        dr = (_r_of(r25, "lo") + _r_of(r25, "hi")) - \
             (_r_of(rb, "lo") + _r_of(rb, "hi"))
        ps = float(r25["Psteel"])
        # convention check: Psteel vs R*I^2 vs 0.5*R*I^2 at the solve drive
        print(f"{cfg}: dRe(Z)_1turn vs no-shell = {dr:.4e} ohm; "
              f"Psteel(FEMM)={ps:.4e} W at 300 A-t "
              f"[R*I^2={dr * 300 ** 2:.4e}, R*I^2/2={dr * 300 ** 2 / 2:.4e}]")
        scale = (a1_pk / 300.0) ** 2
        print(f"   shell eddy loss at ripple: Psteel*(a1_pk/300)^2 = "
              f"{ps * scale * 1e3:.4f} mW (peak-phasor convention) | "
              f"*(a1_rms/300)^2 = {ps * (a1_rms / 300) ** 2 * 1e3:.4f} mW")


# ---------------------------------------------------------------- SIDE part
Z_POS = [-2.5, 0.0, 2.5]
S_VALS = [0.0, S_HOLD, -S_HOLD]
R_WIN = 6.975            # mid of the 0.05 mm coil-shell air gap
Z_ED = [(-13 + k, -12 + k) for k in range(26)]


def side_lua(out_csv: str, fem: str) -> str:
    m = M.magnet
    L: list[str] = []
    w = L.append
    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')
    w('write(handle, "z,s,kind,c1,c2,tot,avg\\n")')

    for z in Z_POS:
        w("newdocument(0)")
        w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
        w('mi_getmaterial("Air")')
        w('mi_getmaterial("Copper")')
        w(f'mi_addmaterial("Magnet", {m.mu_r}, {m.mu_r}, {m.hc:.1f}, '
          "0, 0.667, 0, 0, 1, 0, 0, 0)")
        w('mi_getmaterial("1018 Steel")')
        g = Geo(w)
        g.rect(0, z - m.length / 2, m.radius, z + m.length / 2)
        g.label(m.radius / 2, z, "Magnet", magdir=90, group=1, mesh=0.4)
        for c in M.coils:
            g.rect(c.r_in, c.z_bot, c.r_out, c.z_top)
            w(f'mi_addcircprop("{c.name}", {c.ni:.6g}, 1)')
            g.label(c.r_mean_mm, (c.z_bot + c.z_top) / 2, "Copper",
                    circuit=c.name, turns=1, mesh=0.4)
        s = SHELL
        g.rect(s.r_in, s.z_bot, s.r_out, s.z_top)
        g.label((s.r_in + s.r_out) / 2, 0, "1018 Steel", group=0, mesh=0.2)
        # fine WST air annulus (converge_shell FINE settings)
        g.rect(4.05, -9.7, 4.7, 9.7)
        g.label(4.375, 0, "Air", mesh=0.25)
        w(f"mi_makeABC(7, {M.abc_radius:.6g}, 0, 0, 0)")
        g.label(M.abc_radius * 0.5, M.abc_radius * 0.5, "Air")
        w(f'mark("z={z:.6g} built")')
        w(f'mi_saveas("{z_path(fem)}")')
        for sv in S_VALS:
            for c in M.coils:
                w(f'mi_modifycircprop("{c.name}", 1, {c.ni * sv:.6g})')
            w(f'mark("z={z:.6g} s={sv:.6g} analyze")')
            w("mi_analyze()")
            w('mark("solved")')
            w("mi_loadsolution()")
            w("mo_groupselectblock(1)")
            w("fz = mo_blockintegral(19)")
            w("mo_clearblock()")
            w(f'write(handle, {z:.6g}, ",", {sv:.6g}, '
              f'",Fz,0,0,", fz, ",0\\n")')
            for (z1, z2) in Z_ED:
                w("mo_clearcontour()")
                w(f"mo_addcontour({R_WIN:.6g}, {z1:.6g})")
                w(f"mo_addcontour({R_WIN:.6g}, {z2:.6g})")
                w("tot, avg = mo_lineintegral(0)")     # B.n flux
                w("mo_clearcontour()")
                w(f"mo_addcontour({R_WIN:.6g}, {z1:.6g})")
                w(f"mo_addcontour({R_WIN:.6g}, {z2:.6g})")
                w("t2, a2 = mo_lineintegral(5)")       # (B.n)^2
                w(f'write(handle, {z:.6g}, ",", {sv:.6g}, ",Bn,{z1:.6g},'
                  f'{z2:.6g},", tot, ",", avg, "\\n")')
                w(f'write(handle, {z:.6g}, ",", {sv:.6g}, ",Bn2,{z1:.6g},'
                  f'{z2:.6g},", t2, ",", a2, "\\n")')
            w("mo_clearcontour()")
            w("mo_close()")
        w("mi_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def side_run(timeout: int = 3600):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    base = os.path.join(RESULTS_DIR, "st_w10_slitside")
    lua, fem = base + ".lua", base + ".fem"
    if os.path.exists(SIDE_CSV):
        os.remove(SIDE_CSV)
    with open(lua, "w") as fh:
        fh.write(side_lua(SIDE_CSV, fem))
    subprocess.run(["femm-lua", lua], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    print(f"done: {SIDE_CSV}")


def side_analyze():
    mm = 1e-3
    with open(SIDE_CSV) as fh:
        rows = [dict(kind=r["kind"], z=float(r["z"]), s=float(r["s"]),
                     c1=float(r["c1"]), c2=float(r["c2"]),
                     tot=float(r["tot"]), avg=float(r["avg"]))
                for r in csv.DictReader(fh)]
    kerfs = (0.2, 0.4)
    print("== st_w10 slit side-load (Maxwell tension of the missing wedge, "
          "worst case) ==")
    print("windows: 26 x 1 mm at r=6.975 (shell inner face); "
          "F(w) = sum avgBn2_i*dz*w/(2 mu0)")
    results = {}
    for z in Z_POS:
        for sv in S_VALS:
            bn2 = [(r["c1"], r["c2"], r["avg"]) for r in rows
                   if r["kind"] == "Bn2" and r["z"] == z and r["s"] == sv]
            bn = {(r["c1"]): r["avg"] for r in rows
                  if r["kind"] == "Bn" and r["z"] == z and r["s"] == sv}
            if not bn2:
                continue
            intg = sum(a * (z2 - z1) * mm for (z1, z2, a) in bn2)  # T^2 m
            f_per_w = intg / (2 * MU0)          # N per m of kerf width
            results[(z, sv)] = f_per_w
            pk = max((a, z1) for (z1, _, a) in bn2)
            print(f"z={z:+.1f} s={sv:+.4f}: int Bn^2 dz = {intg:.4e} T^2 m; "
                  f"peak sqrt(avgBn2) = {math.sqrt(pk[0]):.3f} T at "
                  f"z=[{pk[1]:+.0f},{pk[1] + 1:+.0f}]; "
                  + "  ".join(f"F({w}mm)={f_per_w * w * mm * 1e3:.1f}mN"
                              for w in kerfs))
    print("\nprofile Bn(z) at worst case (z=+2.5, s=+S_HOLD) [T]:")
    prof = [(r["c1"], r["avg"]) for r in rows
            if r["kind"] == "Bn" and r["z"] == 2.5 and r["s"] == S_HOLD]
    print("  " + " ".join(f"{z1:+.0f}:{v:+.3f}" for (z1, v) in sorted(prof)))
    worst = max(results.values())
    span = max(results.values()) - min(results.values())
    print(f"\nworst F/w = {worst:.1f} N/m; modulation over stroke/drive = "
          f"{span:.1f} N/m")
    for w in kerfs:
        fw = worst * w * mm
        print(f"kerf {w} mm: single-slit side load {fw * 1e3:.1f} mN "
              f"-> rail friction at mu=0.15: {fw * 0.15 * 1e3:.2f} mN "
              f"(feel threshold 10 mN, cogging 161 mN)")
        dw = 0.1  # +/-0.05 mm kerf tolerance, worst-case mismatch
        ang = math.radians(4.0)  # +/-2 deg placement on both slits
        res = fw * (dw / w) + fw * ang
        print(f"  two-slit residual (kerf mismatch {dw}mm + 4deg): "
              f"{res * 1e3:.2f} mN -> friction {res * 0.15 * 1e3:.3f} mN")
    circ = 2 * math.pi * 7.5
    for nslit, w in ((1, 0.2), (1, 0.4), (2, 0.2), (2, 0.4)):
        loss = nslit * w / circ
        print(f"{nslit} slit(s) x {w} mm kerf: return-path area -{loss * 100:.2f}%"
              f" -> wall bias 1.21 T -> {1.21 / (1 - loss):.3f} T")


# --------------------------------------------------------------------- cli
def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    mode = argv[0]
    if mode == "ac-all":
        ac_all()
    elif mode == "ac-one":
        ac_run_one(argv[1], float(argv[2]))
    elif mode == "ac-analyze":
        ac_analyze()
    elif mode == "side":
        side_run()
    elif mode == "side-analyze":
        side_analyze()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1:])
