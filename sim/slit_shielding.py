"""Transverse-shielding cost of axial slit(s) in the st_w10 shell -- PLANAR 2D.

Why planar: a full-length axial slit breaks axisymmetry, so the r-z solver
cannot represent it. For a *transverse* (x-y) applied field on a long tube the
2D planar cross-section IS the correct idealization (infinite tube; the real
26 mm tube's end effects reduce true shielding somewhat, uniformly across
cases -- we use the planar model for RATIOS against the closed ring and rescale
by the measured axisymmetric-formula S_in ~ 23).

Model
-----
* mi_probdef(0, "millimeters", "planar", 1e-8, DEPTH=26, 30); geometry is the
  shell cross-section: annulus 7.0 <= r <= 8.0 in the x-y plane.
* Uniform applied field B0 = 2 mT (neighbor-field scale at 19.05 mm) imposed
  by a Prescribed-A boundary (BdryFormat 0) on a large square box (half-width
  60 mm): A_z = B0*(y*cos th - x*sin th) -> B = B0*(cos th, sin th). A linear
  in (x,y) is exact on linear elements, so the empty-domain field is uniform
  to mesh precision. The length-unit convention of (x,y) inside FEMM's
  Prescribed-A formula is resolved EMPIRICALLY by `validate` (both candidate
  conventions solved; the one matching 2 mT within 2% is used by `run`).
* LINEAR steel mu_r in {100, 325, 1000}. 325 is the physical case: the wall
  is biased to |B| ~ 1.21 T by the mover's own axial flux, so mT-level
  neighbor fields see the differential mu ~ 325. Do NOT use nonlinear 1018
  here: at a mT-level applied field the solver would sit at initial
  permeability and misrepresent the biased wall. 100/1000 = sensitivity.
* Probes: 0.5 mm square air blocks at origin, (+/-4.75, 0), (0, +/-4.75)
  [Hall-plane radius] and (4.0, 0); local avg B = mo_blockintegral(8|9)/vol.
  Types verified against .femm/app/bin/manual.pdf (mo_blockintegral table:
  8 = integral of Bx, 9 = integral of By, 10 = volume) and cross-checked on
  the empty-domain run where B is known.
* Slits: radial-faced wedges at the given azimuth; angular half-width chosen
  so the gap equals the kerf at mid-wall r = 7.5 (7% narrower at ID, wider at
  OD -- irrelevant at this fidelity). Two-slit cases (azimuth 0 and 180) make
  two independent half-shells.

Field orientations th = angle of applied B from +x (slit at azimuth 0):
  th=0  -> B along the slit axis (slit sits at the natural circumferential
           flux NULL of the shunt path)
  th=90 -> B perpendicular to the slit axis (slit cuts the shunt path at its
           flux MAXIMUM; for two slits BOTH half-ring paths are broken)
  th=45 -> intermediate = what both grid axes see if the slit azimuth is
           rotated 45 deg from the keyboard grid.

Run inside `nix develop`, from the repo root (femm-lua resolves .femm state
from $PWD); subprocess timeouts are set here:
    python sim/slit_shielding.py validate   # 6 small solves, ~1 min
    python sim/slit_shielding.py run        # 41 planar solves
    python sim/slit_shielding.py analyze
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from femm import PROJECT_ROOT, RESULTS_DIR, z_path  # noqa: E402

BASE = os.path.join(RESULTS_DIR, "slit_shield")
VAL_BASE = os.path.join(RESULTS_DIR, "slit_shield_val")

B0 = 2e-3            # T, applied transverse field (neighbor scale at 19 mm)
R_IN, R_OUT = 7.0, 8.0
R_MID = 0.5 * (R_IN + R_OUT)
DEPTH = 26.0         # mm, shell length (z +/-13)
BOX = 60.0           # mm, half-width of prescribed-A square boundary
MUS = (100.0, 325.0, 1000.0)
MU_PHYS = 325.0      # differential mu at the 1.21 T wall bias (measured)
S_IN_MEAS = 23.0     # axisymmetric-formula incoming shielding, closed tube
F_UNSHIELDED_MN = 25.0   # mN worst neighbor force, unshielded victim (DESIGN)
THETAS = (0.0, 45.0, 90.0)

# case -> list of (azimuth_deg, kerf_mm)
CASES = {
    "closed":  [],
    "s1k02":   [(0.0, 0.2)],
    "s1k04":   [(0.0, 0.4)],
    "s2k02":   [(0.0, 0.2), (180.0, 0.2)],
    "s2k04":   [(0.0, 0.4), (180.0, 0.4)],
}

# probe name -> center (mm); 0.5 mm squares
PROBE_H = 0.25
PROBES = {
    "o":   (0.0, 0.0),
    "xp":  (4.75, 0.0),
    "xm":  (-4.75, 0.0),
    "yp":  (0.0, 4.75),
    "ym":  (0.0, -4.75),
    "x4":  (4.0, 0.0),
}

# Prescribed-A (x,y) unit conventions to try: coefficient = B0 * u
CONVS = {"xy_mm": 1e-3, "xy_m": 1.0}


# ------------------------------------------------------------- lua emission
def _header(L, out_csv, cols):
    w = L.append
    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')
    w(f'write(handle, "{cols}\\n")')


def _label(w, x, y, mat, mesh=0.0):
    automesh = 1 if mesh <= 0 else 0
    w(f"mi_addblocklabel({x:.8g}, {y:.8g})")
    w(f"mi_selectlabel({x:.8g}, {y:.8g})")
    w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, "", 0, 0, 0)')
    w("mi_clearselected()")


def _seg(w, x1, y1, x2, y2):
    w(f"mi_addnode({x1:.8g}, {y1:.8g})")
    w(f"mi_addnode({x2:.8g}, {y2:.8g})")
    w(f"mi_addsegment({x1:.8g}, {y1:.8g}, {x2:.8g}, {y2:.8g})")


def _rect(w, x1, y1, x2, y2):
    for (a, b, c, d) in [(x1, y1, x2, y1), (x2, y1, x2, y2),
                         (x2, y2, x1, y2), (x1, y2, x1, y1)]:
        _seg(w, a, b, c, d)


def _ang_diff(a, b):
    d = (a - b) % 360.0
    return min(d, 360.0 - d)


def _ring_angles(slits):
    """Sorted node angles on both circles + list of slit edge angles."""
    edges = []          # per slit: (lo_edge, hi_edge) degrees
    for (phi, kerf) in slits:
        alpha = math.degrees(math.asin(0.5 * kerf / R_MID))
        edges.append(((phi - alpha) % 360.0, (phi + alpha) % 360.0))
    angs = set()
    for (e1, e2) in edges:
        angs.update((round(e1, 6), round(e2, 6)))
    for a in range(0, 360, 45):
        if not any(_ang_diff(a, phi) < math.degrees(
                math.asin(0.5 * k / R_MID)) + 0.5 for (phi, k) in slits):
            angs.add(float(a))
    return sorted(angs), edges


def _emit_ring(w, slits):
    angs, edges = _ring_angles(slits)
    xy = lambda r, a: (r * math.cos(math.radians(a)),
                       r * math.sin(math.radians(a)))
    for a in angs:
        for r in (R_IN, R_OUT):
            x, y = xy(r, a)
            w(f"mi_addnode({x:.8g}, {y:.8g})")
    n = len(angs)
    for i in range(n):
        a1, a2 = angs[i], angs[(i + 1) % n]
        sweep = (a2 - a1) % 360.0
        for r in (R_IN, R_OUT):
            x1, y1 = xy(r, a1)
            x2, y2 = xy(r, a2)
            # maxseg arg = max degrees per element on the arc
            w(f"mi_addarc({x1:.8g}, {y1:.8g}, {x2:.8g}, {y2:.8g}, "
              f"{sweep:.8g}, 1)")
    for (e1, e2) in edges:
        for a in (e1, e2):
            x1, y1 = xy(R_IN, a)
            x2, y2 = xy(R_OUT, a)
            w(f"mi_addsegment({x1:.8g}, {y1:.8g}, {x2:.8g}, {y2:.8g})")
    # steel labels: one per connected sector (between consecutive slits)
    if not slits:
        x, y = xy(R_MID, 90.0)
        _label(w, x, y, "LinSteel", mesh=0.15)
    else:
        phis = sorted(phi for (phi, _) in slits)
        alphas = {phi: math.degrees(math.asin(0.5 * k / R_MID))
                  for (phi, k) in slits}
        for i in range(len(phis)):
            p1, p2 = phis[i], phis[(i + 1) % len(phis)]
            start = p1 + alphas[p1]
            span = (p2 - alphas[p2] - start) % 360.0
            x, y = xy(R_MID, (start + span / 2.0) % 360.0)
            _label(w, x, y, "LinSteel", mesh=0.15)
        for (phi, k) in slits:
            x, y = xy(R_MID, phi)
            _label(w, x, y, "Air", mesh=0.05)


def _emit_solve(w, case, slits, mu, theta, box, u, out_row_prefix, fem):
    """One planar solve: geometry, boundary, probes, post-processing rows."""
    w("newdocument(0)")
    w(f'mi_probdef(0, "millimeters", "planar", 1e-8, {DEPTH:.6g}, 30)')
    w('mi_getmaterial("Air")')
    if slits is not None:
        w(f'mi_addmaterial("LinSteel", {mu:.6g}, {mu:.6g}, 0, 0, 0, 0, '
          "0, 1, 0, 0, 0)")
    # uniform-field boundary: A = A1*x + A2*y (BdryFormat 0)
    th = math.radians(theta)
    a1 = -B0 * math.sin(th) * u
    a2 = B0 * math.cos(th) * u
    w(f'mi_addboundprop("uni", 0, {a1:.9g}, {a2:.9g}, 0, 0, 0, 0, 0, 0)')
    _rect(w, -box, -box, box, box)
    for (sx, sy) in ((0, -box), (box, 0), (0, box), (-box, 0)):
        w(f"mi_selectsegment({sx:.8g}, {sy:.8g})")
    w('mi_setsegmentprop("uni", 0, 1, 0, 0)')
    w("mi_clearselected()")

    if slits is not None:
        _emit_ring(w, slits)
        # mesh-grading circle r=12 between ring and far field
        for (x, y) in ((12, 0), (0, 12), (-12, 0), (0, -12)):
            w(f"mi_addnode({x:.8g}, {y:.8g})")
        for (x1, y1, x2, y2) in ((12, 0, 0, 12), (0, 12, -12, 0),
                                 (-12, 0, 0, -12), (0, -12, 12, 0)):
            w(f"mi_addarc({x1:.8g}, {y1:.8g}, {x2:.8g}, {y2:.8g}, 90, 5)")
        _label(w, 0, 10, "Air", mesh=0.8)       # annulus 8..12
        _label(w, -2.5, 0, "Air", mesh=0.3)     # bore
        _label(w, 0, box / 2.0, "Air", mesh=4.0)  # far field
    else:
        _label(w, 0, box / 2.0, "Air", mesh=2.0)  # empty domain: one region

    for (pname, (px, py)) in PROBES.items():
        _rect(w, px - PROBE_H, py - PROBE_H, px + PROBE_H, py + PROBE_H)
        _label(w, px, py, "Air", mesh=0.1)

    w(f'mark("{case} mu={mu:g} th={theta:g} box={box:g} built")')
    w(f'mi_saveas("{z_path(fem)}")')
    w(f'mark("{case} mu={mu:g} th={theta:g} analyze")')
    w("mi_analyze()")
    w('mark("solved")')
    w("mi_loadsolution()")
    for (pname, (px, py)) in PROBES.items():
        w(f"mo_selectblock({px:.8g}, {py:.8g})")
        w("ibx = mo_blockintegral(8)")
        w("iby = mo_blockintegral(9)")
        w("vol = mo_blockintegral(10)")
        w("mo_clearblock()")
        w(f'write(handle, "{out_row_prefix},{pname},", '
          'ibx, ",", iby, ",", vol, "\\n")')
    w("mo_close()")
    w("mi_close()")


def gen_validate_lua(out_csv, fem):
    L = []
    _header(L, out_csv, "conv,theta,probe,ibx,iby,vol")
    for (conv, u) in CONVS.items():
        for theta in THETAS:
            _emit_solve(L.append, f"empty_{conv}", None, 0.0, theta, BOX, u,
                        f"{conv},{theta:g}", fem)
    L.append("closefile(handle)")
    L.append("quit()")
    return "\n".join(L) + "\n"


def _jobs():
    """(case, mu, theta, box) solve matrix."""
    jobs = []
    for mu in MUS:
        for (case, slits) in CASES.items():
            thetas = (0.0,) if case == "closed" else THETAS
            for th in thetas:
                jobs.append((case, mu, th, BOX))
    jobs.append(("closed", MU_PHYS, 90.0, BOX))    # isotropy check
    jobs.append(("closed", MU_PHYS, 0.0, 2 * BOX))  # boundary doubling check
    return jobs


def gen_run_lua(out_csv, fem, u):
    L = []
    _header(L, out_csv, "case,mu,theta,box,probe,ibx,iby,vol")
    for (case, mu, th, box) in _jobs():
        _emit_solve(L.append, case, CASES[case], mu, th, box, u,
                    f"{case},{mu:g},{th:g},{box:g}", fem)
    L.append("closefile(handle)")
    L.append("quit()")
    return "\n".join(L) + "\n"


def run_femm(lua_text, base, timeout):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    lua, out_csv = base + ".lua", base + ".csv"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua, "w") as f:
        f.write(lua_text)
    subprocess.run(["femm-lua", lua], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    return out_csv


# --------------------------------------------------------------- validation
def read_validate():
    path = VAL_BASE + ".csv"
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({"conv": r["conv"], "theta": float(r["theta"]),
                         "probe": r["probe"], "ibx": float(r["ibx"]),
                         "iby": float(r["iby"]), "vol": float(r["vol"])})
    return rows


def validation_report(verbose=True):
    """Return (winning conv name, u, worst rel err) or raise."""
    rows = read_validate()
    if rows is None:
        raise SystemExit("run `validate` first (no slit_shield_val.csv)")
    best = None
    for (conv, u) in CONVS.items():
        worst = 0.0
        for r in rows:
            if r["conv"] != conv:
                continue
            bx, by = r["ibx"] / r["vol"], r["iby"] / r["vol"]
            th = math.radians(r["theta"])
            ex, ey = B0 * math.cos(th), B0 * math.sin(th)
            err = math.hypot(bx - ex, by - ey) / B0
            worst = max(worst, err)
            if verbose:
                print(f"  {conv} th={r['theta']:3g} {r['probe']:>3}: "
                      f"B=({bx*1e3:+.4f},{by*1e3:+.4f}) mT  "
                      f"want ({ex*1e3:+.4f},{ey*1e3:+.4f})  err {err*100:.2f}%")
        if verbose:
            print(f"  {conv}: worst probe error {worst*100:.2f}% of B0")
        if worst < 0.02 and (best is None or worst < best[2]):
            best = (conv, CONVS[conv], worst)
    if best is None:
        raise SystemExit("NO Prescribed-A convention matched 2 mT within 2% "
                         "-- do not trust the run mode; investigate.")
    print(f"validated: Prescribed-A (x,y) convention '{best[0]}' "
          f"(coeff = B0*{best[1]:g}), worst empty-domain error "
          f"{best[2]*100:.2f}% (< 2% gate)")
    return best


# ----------------------------------------------------------------- analysis
def s_exact_2d(mu, a=R_IN, b=R_OUT):
    """Exact 2D infinite-tube transverse shielding factor."""
    return ((mu + 1) ** 2 * b * b - (mu - 1) ** 2 * a * a) / (4 * mu * b * b)


def read_run():
    path = BASE + ".csv"
    if not os.path.exists(path):
        raise SystemExit("run `run` first (no slit_shield.csv)")
    out = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r["case"], float(r["mu"]), float(r["theta"]),
                   float(r["box"]))
            bx = float(r["ibx"]) / float(r["vol"])
            by = float(r["iby"]) / float(r["vol"])
            out.setdefault(key, {})[r["probe"]] = (bx, by)
    return out


def analyze():
    data = read_run()

    def S_of(key, probe="o"):
        bx, by = data[key][probe]
        return B0 / math.hypot(bx, by)

    print("\n===== planar transverse shielding, st_w10 shell cross-section "
          f"(B0 = {B0*1e3:g} mT) =====")
    print("case      mu   th   S_center  resid_uT  resid_dir  S_worstprobe  "
          "bore spread")
    summary = []
    for key in sorted(data, key=lambda k: (k[0], k[1], k[2], k[3])):
        (case, mu, th, box) = key
        bx, by = data[key]["o"]
        bres = math.hypot(bx, by)
        ang = math.degrees(math.atan2(by, bx))
        svals = {p: B0 / math.hypot(*data[key][p]) for p in data[key]}
        sworst = min(svals.values())
        pw = min(svals, key=svals.get)
        note = " [box2]" if box > BOX else ""
        print(f"{case:9s} {mu:5g} {th:4g}  {B0/bres:8.2f}  {bres*1e6:8.2f} "
              f" {ang:+7.1f}   {sworst:8.2f} ({pw})   "
              f"{min(svals.values())/max(svals.values()):.3f}{note}")
        summary.append({"case": case, "mu": mu, "theta": th, "box": box,
                        "S_center": B0 / bres, "resid_uT": bres * 1e6,
                        "resid_dir_deg": ang, "S_worst_probe": sworst,
                        "worst_probe": pw})

    # ---- model consistency checks
    print("\n---- consistency checks ----")
    for mu in MUS:
        k = ("closed", mu, 0.0, BOX)
        if k in data:
            print(f"closed mu={mu:g}: planar S = {S_of(k):.2f} vs exact-2D "
                  f"formula {s_exact_2d(mu):.2f} vs thin-wall mu*t/2R+1 = "
                  f"{mu / (2 * R_MID) + 1:.2f}")
    k0, k90 = ("closed", MU_PHYS, 0.0, BOX), ("closed", MU_PHYS, 90.0, BOX)
    if k90 in data:
        print(f"closed isotropy: S(0)={S_of(k0):.2f} S(90)={S_of(k90):.2f} "
              f"(ratio {S_of(k0)/S_of(k90):.3f})")
    kb2 = ("closed", MU_PHYS, 0.0, 2 * BOX)
    if kb2 in data:
        print(f"boundary doubling: S(box60)={S_of(k0):.2f} "
              f"S(box120)={S_of(kb2):.2f} (ratio {S_of(k0)/S_of(kb2):.3f})")
    print(f"axisym-formula S_in ~ {S_IN_MEAS:g} vs planar closed "
          f"S={S_of(('closed', MU_PHYS, 0.0, BOX)):.2f} -> "
          f"rescale factor {S_IN_MEAS/S_of(('closed', MU_PHYS, 0.0, BOX)):.3f}"
          " (applied below)")

    # ---- keyboard translation at the physical mu
    s_closed = S_of(("closed", MU_PHYS, 0.0, BOX))
    scale = S_IN_MEAS / s_closed

    def s_eff(case, th):
        return S_of((case, MU_PHYS, th, BOX)) * scale

    print(f"\n---- keyboard numbers (mu={MU_PHYS:g}, S_eff = planar S x "
          f"{scale:.3f}; victim force = {F_UNSHIELDED_MN:g} mN / S_eff; "
          "feel threshold 10 mN) ----")
    print("kerf removed from DC force-path circumference: "
          + ", ".join(f"{c}: {sum(k for _, k in s) / (2*math.pi*R_MID)*100:.2f}%"
                      for c, s in CASES.items() if s))
    hdr = (f"{'case':9s} {'orient':26s} {'S_eff':>6s} {'F_nb[mN]':>9s} "
           f"{'%feel':>6s} {'%0.8N':>6s}")
    print(hdr)

    def row(case, label, seff):
        f = F_UNSHIELDED_MN / seff
        print(f"{case:9s} {label:26s} {seff:6.1f} {f:9.2f} "
              f"{f/10*100:6.0f} {f/800*100:6.2f}")

    row("closed", "any", s_eff("closed", 0.0))
    for case in ("s1k02", "s1k04", "s2k02", "s2k04"):
        for (th, lab) in ((0.0, "B || slit axis (best)"),
                          (45.0, "B at 45 deg"),
                          (90.0, "B perp slit axis (worst)")):
            row(case, lab, s_eff(case, th))
        # grid choices: victim has orthogonal neighbors on BOTH axes
        s_grid = min(s_eff(case, 0.0), s_eff(case, 90.0))
        s_diag = s_eff(case, 45.0)
        row(case, "slit on grid axis: binding", s_grid)
        row(case, "slit at 45 to grid: both", s_diag)
    print("\nnotes:")
    print(" * grid-aligned slit: one neighbor pair sees th=0, the other "
          "th=90 -> the th=90 number binds.")
    print(" * slit at 45 deg to grid: all 4 orthogonal neighbors see th=45; "
          "the th=90 worst case moves to the 4 diagonal neighbors at "
          "26.9 mm, whose unshielded drive is ~0.25-0.35x of 25 mN "
          "(dipole-scaling estimate ~6-9 mN).")
    print(" * source leakage (x18) is unaffected by the slit to first "
          "order; only the victim-response 1/S_in part changes.")
    print(" * sensor: the dominant 0.18 mT neighbor field at the Hall is "
          "AXIAL (unshielded by the tube already) -> slit-independent. The "
          "transverse component (same order) rises from ~S_eff-shielded to "
          "worst-case slit levels; for an axial-sensing Hall this enters "
          "only via cross-axis sensitivity (~1-2%) -> stays um-class.")

    path = BASE + "_summary.csv"
    with open(path, "w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(summary)
    print(f"\nwrote {path}")


# ---------------------------------------------------------------------- cli
def main(argv):
    if not argv or argv[0] == "analyze":
        analyze()
        return
    if argv[0] == "validate":
        run_femm(gen_validate_lua(VAL_BASE + ".csv", VAL_BASE + ".fem"),
                 VAL_BASE, timeout=900)
        validation_report(verbose=True)
        return
    if argv[0] == "run":
        conv, u, err = validation_report(verbose=False)
        run_femm(gen_run_lua(BASE + ".csv", BASE + ".fem", u),
                 BASE, timeout=3000)
        print(f"done: {BASE}.csv")
        return
    raise SystemExit(f"unknown mode {argv}")


if __name__ == "__main__":
    main(sys.argv[1:])
