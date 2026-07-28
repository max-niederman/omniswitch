"""Convergence qualification of the thick-shell single-magnet variants
st_w05L (0.5 mm wall, long coil) and st_w10 (1.0 mm wall) as THE design.

These are the two surviving shell thicknesses from sim/explore_steel_shell.py
(pot/washer variants rejected for ~0.8 N cogging; st_w03's 0.3 mm shell is
crosstalk-saturated per sim/qualify_crosstalk.py). This file qualifies them on
every axis that killed or wounded the other finalists:

  (a) fine-mesh stroke sweep    -> hold power, cogging (+/-s differencing AND
                                   s=0 rows, z-antisymmetry test), gain
                                   flatness/saturation, tau, package budget
  (b) crosstalk source leakage  -> contour-window shell/bare ratios at
                                   r=15.05/19.05/23.05 (qualify_crosstalk
                                   method + calibrated analytic propagation),
                                   plus per-slab shell |B| block integrals to
                                   justify an unsaturated-wall mu_r_eff for
                                   incoming shielding S_in
  (c) eddy lag                  -> sim/explore_steel_ac.py harmonic study
                                   (its VARIANTS already include both walls)
  (d) demag                     -> sim/qualify_demag.py sub-block machinery,
                                   st_w10, N45SH/N38UH, T=80/100/120,
                                   s_peak = s_hold_worst / 0.30 from (a)

Run inside `nix develop` from the repo root, ALWAYS under `timeout`:
    python sim/converge_shell.py sweep st_w05L        # ~10 min (77 solves)
    python sim/converge_shell.py sweep st_w10
    python sim/converge_shell.py report st_w05L       # offline re-analysis
    python sim/converge_shell.py xtalk st_w05L        # 1 solve + contours
    python sim/converge_shell.py xtalk st_w10
    python sim/converge_shell.py xtalk-analyze        # offline
    python sim/converge_shell.py ac st_w10            # 30 harmonic solves
    python sim/converge_shell.py ac st_w05L
    python sim/converge_shell.py demag st_w10         # 36 solves (needs sweep)
    python sim/converge_shell.py demag-analyze        # offline, vs st_w03
"""

from __future__ import annotations

import csv
import dataclasses
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from femm import MU0, RESULTS_DIR, Coil, Magnet, Steel, run_sweep  # noqa: E402
from analyze import evaluate, summarize, scale_for_force, total_power, \
    coil_tau, HOLD_DUTY  # noqa: E402
from designs import ShellModel  # noqa: E402
import qualify_crosstalk as qc  # noqa: E402
import qualify_demag as qd  # noqa: E402
from magnets import GRADES  # noqa: E402

TAG = "cvg"
Z_FINE = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
NI_SCALES = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
mm = 1e-3

# ---------------------------------------------------------------- variants
# Geometry copied EXACTLY from sim/explore_steel_shell.py VARIANTS; fine-mesh
# settings copied from its FINAL st_w03_fine (mesh_air 0.4 + 0.25 mm air
# annulus in the bore gap = the WST weighting-gradient region).
FINE = dict(mesh_air=0.4, fine_annulus=(4.05, 4.7, -9.7, 9.7, 0.25))
MAG = Magnet(radius=4.0, length=12.0, z_center=0.0)

VARIANTS: dict[str, ShellModel] = {
    "st_w05L": ShellModel(
        name="st_w05L",
        magnet=MAG,
        coils=[Coil("lo", 4.75, 7.45, -14.5, -0.5, ni=-300.0),
               Coil("hi", 4.75, 7.45, 0.5, 14.5, ni=+300.0)],
        steels=[Steel("shell", 7.5, 8.0, -15.0, 15.0)],
        **FINE),
    "st_w10": ShellModel(
        name="st_w10",
        magnet=MAG,
        coils=[Coil("lo", 4.75, 6.95, -12.5, -0.5, ni=-300.0),
               Coil("hi", 4.75, 6.95, 0.5, 12.5, ni=+300.0)],
        steels=[Steel("shell", 7.0, 8.0, -13.0, 13.0)],
        **FINE),
}

# radial budget (mm): magnet 0..4.0 fixed, OD/2 = 8.0 fixed
BUDGET = {
    "st_w03":  dict(gap=(4.0, 4.75), coil=(4.75, 7.65), shell=(7.7, 8.0)),
    "st_w05L": dict(gap=(4.0, 4.75), coil=(4.75, 7.45), shell=(7.5, 8.0)),
    "st_w10":  dict(gap=(4.0, 4.75), coil=(4.75, 6.95), shell=(7.0, 8.0)),
}


def package_len(m: ShellModel, stroke_half=2.5):
    lo = min(min(c.z_bot for c in m.coils), min(s.z_bot for s in m.steels),
             -stroke_half - m.magnet.length / 2)
    hi = max(max(c.z_top for c in m.coils), max(s.z_top for s in m.steels),
             stroke_half + m.magnet.length / 2)
    return hi - lo + 2.0     # + 2 mm end caps


def sweep_csv(name: str) -> str:
    return os.path.join(RESULTS_DIR, f"{name}_{TAG}.csv")


def load_rows(name: str):
    with open(sweep_csv(name)) as f:
        return [{k: float(v) for k, v in row.items()}
                for row in csv.DictReader(f)]


# ------------------------------------------------------------- (a) sweep
def do_sweep(name: str, timeout: int = 3600):
    m = VARIANTS[name]
    rows = run_sweep(m, Z_FINE, NI_SCALES, tag=TAG, timeout=timeout)
    report(name, rows)


def report(name: str, rows=None):
    m = VARIANTS[name]
    rows = rows or load_rows(name)
    ev = evaluate(m, rows)
    print(summarize(ev))
    p_worst = ev["power_worst_W"]
    print(f"fill 0.6 worst hold power: {p_worst:.3f} W ; "
          f"fill 0.5 (honest): {p_worst * 1.2:.3f} W")
    s_hold = max(max(abs(v) for v in ev["scale_up"].values()),
                 max(abs(v) for v in ev["scale_down"].values()))
    print(f"s_hold_worst = {s_hold:.4f} -> s_peak = s_hold/0.30 = "
          f"{s_hold / 0.30:.4f}")

    byz: dict[float, dict[float, float]] = {}
    for r in rows:
        byz.setdefault(r["z"], {})[r["ni_scale"]] = r["Fz"]

    print("\nforce gain flatness / saturation (antisym +/-s differencing):")
    print("  z[mm]    g05       g1        g2      sat%=1-g2/g05  "
          "even05[N]  even1[N]   s0[N]")
    for z in sorted(byz):
        f = byz[z]
        g05 = (f[0.5] - f[-0.5]) / 1.0
        g1 = (f[1.0] - f[-1.0]) / 2.0
        g2 = (f[2.0] - f[-2.0]) / 4.0
        e05 = (f[0.5] + f[-0.5]) / 2.0     # even part = cogging + even error
        e1 = (f[1.0] + f[-1.0]) / 2.0
        print(f"  {z:5.2f}  {g05:8.4f}  {g1:8.4f}  {g2:8.4f}  "
              f"{100 * (1 - g2 / g05):9.2f}%   {e05:+8.4f}  {e1:+8.4f}  "
              f"{f[0.0]:+8.4f}")

    print("\ncogging z-(anti)symmetry (real cogging must be odd in z):")
    print("  z[mm]   s0(z)      s0(-z)     antisym=(s0(z)-s0(-z))/2   "
          "sym-resid=(s0(z)+s0(-z))/2")
    cogmax = 0.0
    for z in sorted(byz):
        if z <= 0 or -z not in byz:
            continue
        a, b = byz[z][0.0], byz[-z][0.0]
        anti, symr = (a - b) / 2, (a + b) / 2
        cogmax = max(cogmax, abs(anti))
        print(f"  {z:5.2f}  {a:+8.4f}  {b:+8.4f}   {anti:+8.4f}"
              f"                    {symr:+8.4f}")
    print(f"max |antisym cogging| = {cogmax:.4f} N "
          f"({'firmware-compensable' if cogmax < 0.3 else 'TOO BIG'} vs 0.3 N)")

    tau = coil_tau(m, rows)
    print(f"\ntau (flux columns) = {tau * 1e6:.0f} us ; rise to hold at "
          f"{HOLD_DUTY:.0%} duty = {tau * math.log(1 / (1 - HOLD_DUTY)) * 1e6:.0f} us")

    b = BUDGET[name]
    cw = b["coil"][1] - b["coil"][0]
    sw = b["shell"][1] - b["shell"][0]
    print(f"\nradial budget (mm, OD 16.0): magnet 0-4.0 | gap 4.0-"
          f"{b['gap'][1]:.2f} (0.75: bore tube+clearance) | coil "
          f"{b['coil'][0]:.2f}-{b['coil'][1]:.2f} ({cw:.2f}) | air "
          f"{b['coil'][1]:.2f}-{b['shell'][0]:.2f} | shell "
          f"{b['shell'][0]:.2f}-{b['shell'][1]:.2f} ({sw:.2f})")
    plen = package_len(m)
    print(f"package length: {plen:.1f} mm "
          f"({'PASS' if plen <= 35 else 'FAIL'} <= 35)")


# ------------------------------------------------------------ (b) crosstalk
# qc.Design entries (analytic source = the bare single D8x12, same as st_w03;
# s08 refined after the fine sweep exists).
def _qc_design(name: str) -> qc.Design:
    m = VARIANTS[name]
    sh = m.steels[0]
    return qc.Design(
        name,
        magnets=[qc.Mag(4 * mm, 12 * mm, 0.0, +1)],
        coils=[qc.Cl(c.name, c.r_in * mm, c.r_out * mm, c.z_bot * mm,
                     c.z_top * mm, c.ni) for c in m.coils],
        s08=0.35,                        # placeholder; refined from sweep CSV
        sensor_z=-9.5 * mm,
        mover_bot_off=-6.0 * mm,
        shell=(sh.r_in * mm, sh.r_out * mm, sh.z_top * mm),
        note=f"single D8x12 + {sh.r_out - sh.r_in:.1f} mm shell")


for _n in VARIANTS:
    qc.DESIGNS[_n] = _qc_design(_n)
    qc.FEMM_BASE[_n] = f"xtalk_src_{_n}"


def gen_src_slab_lua(d: qc.Design, out_csv: str, fem: str,
                     slab_mm: float = 2.0) -> str:
    """qc.gen_src_lua clone with the shell split into axial slabs (own group
    each) so per-slab avg Br/Bz block integrals report the wall |B|
    (saturation check for mu_r_eff). Emits qc's exact contour windows so
    femm_ratio/k_src machinery applies unchanged."""
    L: list[str] = []
    w = L.append
    qc._lua_header(L, out_csv)
    w("newdocument(0)")
    w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w(f'mi_addmaterial("Magnet", {qc.MUR}, {qc.MUR}, {qc.HC_STR}, '
      "0, 0.667, 0, 0, 1, 0, 0, 0)")
    w('mi_getmaterial("1018 Steel")')
    g = qc._Geo(w)
    for mg in d.magnets:
        r_, l_ = mg.a / mm, mg.L / mm
        g.rect(0, -l_ / 2, r_, l_ / 2)
        g.label(r_ / 2, 0, "Magnet", magdir=90, group=1, mesh=0.4)
    r1, r2, zh = (v / mm for v in d.shell)
    n = int(round(2 * zh / slab_mm))
    zed = [-zh + 2 * zh * k / n for k in range(n + 1)]
    slabs = []
    for k in range(n):
        g.rect(r1, zed[k], r2, zed[k + 1])
        zc = 0.5 * (zed[k] + zed[k + 1])
        g.label((r1 + r2) / 2, zc, "1018 Steel", group=10 + k, mesh=0.2)
        slabs.append((10 + k, zc))
    # far-field mesh control identical to qc.gen_src_lua
    g.rect(14.5, -26, 24, 26)
    g.label(19.25, 25, "Air", mesh=1.0)
    g.rect(0, 14, 2.6, 25)
    g.label(1.3, 20, "Air", mesh=0.5)
    w("mi_makeABC(7, 60, 0, 0, 0)")
    g.label(30, 30, "Air")
    w('mark("built")')
    w(f'mi_saveas("{qc.z_path(fem)}")')
    w('mark("analyze")')
    w("mi_analyze()")
    w('mark("solved")')
    w("mi_loadsolution()")
    w("mo_groupselectblock(1)")
    w("fz = mo_blockintegral(19)")
    w("mo_clearblock()")
    w('write(handle, "Fznoise,0,0,0,", fz, ",0\\n")')
    for gid, zc in slabs:
        w(f"mo_groupselectblock({gid})")
        w("ibr = mo_blockintegral(8)")
        w("ibz = mo_blockintegral(9)")
        w("vol = mo_blockintegral(10)")
        w("mo_clearblock()")
        # kind=shellB, c1=z_center, c2=avgBr, c3=avgBz (T)
        w(f'write(handle, "shellB,{zc:.6g},", ibr/vol, ",", ibz/vol, '
          '",0,0\\n")')
    qc._emit_contours(w, qc.contour_defs(axis_windows=True))
    w("mo_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def do_xtalk(name: str):
    d = qc.DESIGNS[name]
    base = os.path.join(RESULTS_DIR, qc.FEMM_BASE[name])
    qc.run_femm(gen_src_slab_lua(d, base + ".csv", base + ".fem"), base)
    print(f"done: {base}.csv")


# 1018 matlib BH points (B[T], H[A/m]) for differential-mu at the wall bias
BH_1018 = [(0.0, 0.0), (0.2503, 238.732), (0.925, 795.775), (1.25, 1591.55),
           (1.39, 2387.325), (1.525, 3978.875), (1.71, 7957.75),
           (1.87, 15915.5), (1.955, 23873.25), (2.02, 39788.75),
           (2.11, 79577.5), (2.225, 159155.0), (2.43, 318310.0)]


def mur_diff_at(b: float) -> float:
    """Differential relative permeability dB/dH/mu0 of 1018 at |B|=b."""
    for (b1, h1), (b2, h2) in zip(BH_1018, BH_1018[1:]):
        if b <= b2:
            return (b2 - b1) / (h2 - h1) / MU0
    return 1.0


def _ksrc_br(tag: str, field_rz):
    """median |FEMM/analytic| over the Br windows only, per radius + overall."""
    rows = qc.read_femm(tag)
    per_r: dict[float, list[float]] = {}
    for r in rows:
        if r["kind"] != "Br":
            continue
        area = 2 * np.pi * r["c1"] * mm * (r["c3"] - r["c2"]) * mm
        _, an = qc.window_avg(field_rz, "Br", r["c1"], r["c2"], r["c3"])
        if abs(an) > 1e-9:
            per_r.setdefault(r["c1"], []).append(abs(r["tot"] / area) / abs(an))
    med_all = float(np.median([v for lst in per_r.values() for v in lst]))
    return med_all, {rr: float(np.median(v)) for rr, v in per_r.items()}


def _shell_ratio(tag: str, ref: str = "st_w03_bare"):
    """|tot| ratio of Br windows vs the bare-magnet run, per radius."""
    a, b = qc.read_femm(tag), qc.read_femm(ref)
    key = lambda r: (r["kind"], r["c1"], r["c2"], r["c3"])
    bd = {key(r): r for r in b}
    per_r: dict[float, list[float]] = {}
    for r in a:
        if r["kind"] != "Br":
            continue
        rb = bd.get(key(r))
        if rb and abs(rb["tot"]) > 1e-12:
            per_r.setdefault(r["c1"], []).append(abs(r["tot"]) / abs(rb["tot"]))
    return {rr: (float(np.median(v)), min(v), max(v))
            for rr, v in per_r.items()}


def _refine_s08(name: str):
    if os.path.exists(sweep_csv(name)):
        ev = evaluate(VARIANTS[name], load_rows(name))
        s_hold = max(max(abs(v) for v in ev["scale_up"].values()),
                     max(abs(v) for v in ev["scale_down"].values()))
        qc.DESIGNS[name].s08 = s_hold


def xtalk_analyze():
    D = qc.PITCH
    zgrid = np.linspace(-qc.STROKE, qc.STROKE, 5)
    bare = qc.DESIGNS["st_w03"]      # same magnet, analytic source

    for name in ("st_w03", "st_w05L", "st_w10"):
        _refine_s08(name) if name in VARIANTS else None
        d = qc.DESIGNS[name]
        tag = name
        print(f"\n===== {name} ({d.note}) =====")
        rat = _shell_ratio(tag)
        for rr in sorted(rat):
            med, lo, hi = rat[rr]
            print(f"  shell/bare leakage |Br| at r={rr:g}: median {med:.4f} "
                  f"(range {lo:.4f}..{hi:.4f}) -> attenuation x{1 / med:.1f}")
        ksrc, ksrc_r = _ksrc_br(tag, qc.src_field(bare, 0.0))
        print(f"  k_src (|FEMM/analytic-bare| Br windows): median {ksrc:.4f}; "
              f"per radius { {k: round(v, 4) for k, v in ksrc_r.items()} }")

        # shell wall B (slab block integrals)
        rows = qc.read_femm(tag)
        slabs = [(r["c1"], r["c2"], r["c3"]) for r in rows
                 if r["kind"] == "shellB"]
        if slabs:
            bmax = max(math.hypot(br, bz) for (_, br, bz) in slabs)
            mur = mur_diff_at(bmax)
            sh = VARIANTS[name].steels[0] if name in VARIANTS else \
                Steel("shell", 7.7, 8.0, -13, 13)
            t = sh.r_out - sh.r_in
            rmean = 0.5 * (sh.r_in + sh.r_out)
            s_in = 1 + mur * t / (2 * rmean)
            print(f"  shell slab |B|: max {bmax:.3f} T "
                  f"(slab avgs; peak at z={max(slabs, key=lambda s: math.hypot(s[1], s[2]))[0]:+.1f} mm)")
            print("   slab z, avgBr, avgBz [T]: "
                  + "  ".join(f"({z:+.0f}:{math.hypot(br, bz):.2f})"
                              for z, br, bz in slabs))
            print(f"  differential mu_r of 1018 at that bias: {mur:.0f} -> "
                  f"transverse incoming shield S_in ~ 1 + mu_t/(2R) = "
                  f"{s_in:.0f} (infinite-tube formula; axial ~ unshielded)")

        # propagate to neighbor force + sensor with k_src (folds shell)
        f_worst, fx_worst, at = 0.0, 0.0, (0, 0)
        for z_src in zgrid:
            f = qc.scaled(qc.src_field(bare, z_src), ksrc)
            for z_nb in zgrid:
                fx, fz = qc.force_on_neighbor(f, bare, z_nb, D)
                if abs(fz) > f_worst:
                    f_worst, at = abs(fz), (z_src, z_nb)
                fx_worst = max(fx_worst, abs(fx))
        print(f"  neighbor mover force (k_src-scaled charge model): "
              f"|Fz| = {f_worst * 1e3:.1f} mN at z_src={at[0] * 1e3:+.1f}, "
              f"z_nb={at[1] * 1e3:+.1f}; |Fx| up to {fx_worst * 1e3:.1f} mN")
        print(f"    = {f_worst / 0.8 * 100:.1f}% of 0.8 N spec; "
              f"{f_worst / 0.010 * 100:.0f}% of 10 mN feel threshold")

        # sensor corruption (axial Bz at neighbor sensor, unshielded)
        own = qc.src_field(bare, 0.0)
        eps = 0.05 * mm
        _, bz_p = qc.src_field(bare, +eps)(np.array([0.0]),
                                           np.array([d.sensor_z]))
        _, bz_m = qc.src_field(bare, -eps)(np.array([0.0]),
                                           np.array([d.sensor_z]))
        slope = (bz_p[0] - bz_m[0]) / (2 * eps)
        bzx = []
        for z0 in zgrid:
            f = qc.scaled(qc.src_field(bare, z0), ksrc)
            _, _, bz = qc.B_xyz(f, np.array([D]), np.array([0.0]),
                                np.array([d.sensor_z]))
            bzx.append(bz[0])
        bzx = np.array(bzx)
        dbz = bzx.max() - bzx.min()
        fc = qc.scaled(qc.src_coil_field(d, d.s08), ksrc)
        _, _, bz_c = qc.B_xyz(fc, np.array([D]), np.array([0.0]),
                              np.array([d.sensor_z]))
        print(f"  neighbor Hall (axial, unshielded): B_z = "
              f"{bzx[2] * 1e3:+.3f} mT static, p-p over src stroke "
              f"{dbz * 1e6:.0f} uT; own gain {slope * 1e-3:.1f} mT/mm "
              f"-> pos corruption {abs(dbz / slope) * 1e6:.1f} um")
        print(f"    coil-drive xtalk at sensor (src holding 0.8 N, "
              f"s08={d.s08:.3f}): {bz_c[0] * 1e6:+.1f} uT -> "
              f"{abs(2 * bz_c[0] / slope) * 1e6:.2f} um")


# ------------------------------------------------------------------ (c) ac
def do_ac(name: str):
    import explore_steel_ac as ac
    m = ac.VARIANTS[name]           # same geometry (explore_steel_shell)
    rows = ac.run_ac(m)
    ac.analyze_ac(m, rows)
    ac_summary(name)


def ac_summary(name: str):
    path = os.path.join(RESULTS_DIR, f"{name}_stlac.csv")
    with open(path) as f:
        rows = list(csv.DictReader(f))
    print(f"== {name} equivalent first-order lag (hi coil) ==")
    for cfg in ("steel529", "steel100", "slit529"):
        rr = sorted([r for r in rows if r["config"] == cfg],
                    key=lambda r: float(r["freq"]))
        ref = next(r for r in rr if float(r["freq"]) == 0.0)
        lam0, arg0 = float(ref["absLam_hi"]), float(ref["argLam_hi"])
        taus = []
        for r in rr:
            fq = float(r["freq"])
            if fq == 0:
                continue
            phi = float(r["argLam_hi"]) - arg0
            while phi > math.pi:
                phi -= 2 * math.pi
            while phi < -math.pi:
                phi += 2 * math.pi
            taus.append((fq, math.tan(-phi) / (2 * math.pi * fq) * 1e6,
                         float(r["absLam_hi"]) / lam0))
        tline = "  ".join(f"{f:.0f}Hz:{t:.1f}us(|l|={a:.4f})"
                          for f, t, a in taus)
        print(f"  {cfg}: {tline}")


# ---------------------------------------------------------------- (d) demag
def do_demag(name: str = "st_w10"):
    m = VARIANTS[name]
    ev = evaluate(m, load_rows(name))
    s_hold = max(max(abs(v) for v in ev["scale_up"].values()),
                 max(abs(v) for v in ev["scale_down"].values()))
    sh = m.steels[0]
    d = qd.Design(
        name=name,
        parts=[qd.Part("A", 0.0, 12.0, +1, facing=0)],
        coils=list(m.coils),
        shell=(sh.r_in, sh.r_out, sh.z_bot, sh.z_top),
        s_hold=round(s_hold, 4), s_peak=round(s_hold / 0.30, 4))
    qd.DESIGNS[name] = d
    print(f"{name}: s_hold={d.s_hold} s_peak={d.s_peak} "
          f"(peak NI/coil = {300 * d.s_peak:.0f} A-t)")
    qd.run_design(d, grades=["N45SH", "N38UH"], temps=[80.0, 100.0, 120.0],
                  zs=[0.0, 2.5], tag="demag")
    print(f"done {name}")


def demag_analyze():
    for name in ("st_w03", "st_w10"):
        path = os.path.join(RESULTS_DIR, f"{name}_demag_blocks.csv")
        if not os.path.exists(path):
            print(f"{name}: no demag data")
            continue
        rows = qd.margin_rows(path, name)
        mins: dict = {}
        mins0: dict = {}
        for r in rows:
            k = (r["grade"], r["temp"])
            if k not in mins or r["margin_peak"] < mins[k]["margin_peak"]:
                mins[k] = r
            if r["s"] == 0.0 and (k not in mins0
                                  or r["margin_peak"] < mins0[k]["margin_peak"]):
                mins0[k] = r
        print(f"\n== {name} min demag margin (1.5x peaking) ==")
        for (g, T) in sorted(mins):
            r = mins[(g, T)]
            r0 = mins0[(g, T)]
            print(f"  {g:6s} T={T:5.0f}: margin_peak={r['margin_peak']:6.2f} "
                  f"(avg {r['margin_avg']:6.2f}; s=0 only "
                  f"{r0['margin_peak']:6.2f}) worst {r['part']}/{r['slab']}/"
                  f"{r['zone']} z={r['z']} s={r['s']:+.3f} "
                  f"avgBz={r['avg_bz']:+.3f} T")

    # hold power per grade/temp for st_w10 from its demag force gains
    fpath = os.path.join(RESULTS_DIR, "st_w10_demag_force.csv")
    if os.path.exists(fpath) and "st_w10" in qd.DESIGNS:
        d = qd.DESIGNS["st_w10"]
    else:
        d = None
    if d is None and os.path.exists(fpath):
        # reconstruct s_peak from sweep
        ev = evaluate(VARIANTS["st_w10"], load_rows("st_w10"))
        s_hold = max(max(abs(v) for v in ev["scale_up"].values()),
                     max(abs(v) for v in ev["scale_down"].values()))
        d = qd.Design(name="st_w10",
                      parts=[qd.Part("A", 0.0, 12.0, +1, facing=0)],
                      coils=list(VARIANTS["st_w10"].coils),
                      s_hold=round(s_hold, 4),
                      s_peak=round(s_hold / 0.30, 4))
    if d is not None and os.path.exists(fpath):
        force = {}
        with open(fpath) as f:
            for r in csv.DictReader(f):
                force[(r["grade"], float(r["temp"]), float(r["z"]),
                       float(r["s"]))] = float(r["Fz"])
        print("\n== st_w10 hold power per grade/temp (0.8 N worst dir, "
              "fill 0.6 | fill 0.5) ==")
        for g in ("N45SH", "N38UH"):
            cells = []
            for T in (80.0, 100.0, 120.0):
                gains = []
                for z in (0.0, 2.5):
                    fp = force[(g, T, z, d.s_peak)]
                    fm = force[(g, T, z, -d.s_peak)]
                    gains.append((fp - fm) / (2 * d.s_peak))
                s_hold_t = 0.8 / min(gains)
                p06 = sum(c.power(c.ni * s_hold_t) for c in d.coils)
                cells.append(f"T={T:.0f}: {p06:.3f}|{p06 * 1.2:.3f} W")
            print(f"  {g}: " + "   ".join(cells))


# ------------------------------------------------------------------ verdict
# Qualification results (2026-07-28 runs; CSVs in results/*_cvg.csv,
# results/xtalk_src_st_w{05L,10}.csv, results/st_w{05L,10}_stlac.csv,
# results/st_w10_demag_*.csv):
#
#                              st_w05L        st_w10        (st_w03 ref)
#   hold worst, N52-sim Br1.43 1.00/1.20 W    1.27/1.53 W   1.01/1.22 W  (fill .6/.5)
#   cogging (z-antisym, real)  0.039 N        0.161 N       0.068 N     all < 0.3 OK
#   sym-resid (mesh noise)     0.033 N        0.014 N       -
#   F-I saturation |1-g2/g05|  0.01%          0.00%         <2%
#   gain center/end [N/unit]   2.377/2.186    2.818/2.535   2.526/2.289
#   tau                        170 us         137 us        162 us
#   package (len x OD)         32.0 x 16.0    28.0 x 16.0   28.0 x 16.0
#   leak vs bare @19.05 mm     0.229 (x4.4)   0.055 (x18)   0.486 (x2.1)
#   neighbor force (worst)     98 mN          25 mN         204 mN
#   neighbor Hall corruption   13.3 um        3.3 um        27.7 um
#   shell peak |B| (slab avg)  1.85 T SAT     1.21 T UNSAT  (>2 T SAT)
#   S_in transverse (mu_diff)  ~2 (mur 16)    ~23 (mur 325) ~2
#   eddy lag equiv             3.2 us         5.7 us        1.9 us  (~prop t, not t^2)
#   demag N45SH/N38UH @100C    not run        1.40/1.94     1.38/1.91
#
# VERDICT: st_w10. Crosstalk was the killing axis and only the 1.0 mm wall
# is genuinely unsaturated -> 18x source attenuation (25 mN residual, 3% of
# spec) + real incoming shielding; costs +27% hold power vs st_w05L and a
# 0.16 N (compensable, centering) cogging; demag unchanged vs st_w03
# (N45SH passes 1.3 @100C, N38UH for 120 C headroom). st_w05L's 0.5 mm wall
# saturates at 1.85 T (x4.4 only, 98 mN ~ 10x feel threshold) — rejected.
# Do NOT thin below ~1.0 mm: wall bias 1.21 T already; 0.8 mm would sit
# ~1.5 T and start saturating (cold N52-class magnet is the worst case).
# Grade hold power (worst dir incl. cogging, Br-scaled): N45SH 20C
# 1.42/1.70 W, 100C 1.63/1.95 W; N38UH 100C 1.81/2.17 W (fill .6/.5).


# -------------------------------------------------------------------- cli
def main(argv):
    if not argv:
        raise SystemExit(__doc__)
    mode = argv[0]
    if mode == "sweep":
        do_sweep(argv[1], timeout=int(argv[2]) if len(argv) > 2 else 3600)
    elif mode == "report":
        report(argv[1])
    elif mode == "xtalk":
        do_xtalk(argv[1])
    elif mode == "xtalk-analyze":
        xtalk_analyze()
    elif mode == "ac":
        do_ac(argv[1])
    elif mode == "ac-summary":
        ac_summary(argv[1])
    elif mode == "demag":
        do_demag(argv[1] if len(argv) > 1 else "st_w10")
    elif mode == "demag-analyze":
        demag_analyze()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main(sys.argv[1:])
