"""Inter-switch magnetic crosstalk qualification at 19.05 mm key pitch.

FEMM is axisymmetric and cannot model two parallel actuators, so the two-body
quantities (field at the neighbor's sensor, force injected on the neighbor's
mover, coil-coil transformer coupling) are computed ANALYTICALLY here, with an
exact numeric field model (no far-field multipole truncation -- at 19 mm the
movers are NOT in their far field):

  * axially magnetized cylinder == solenoid sheet K = Hc = Br/(mu0*mu_r),
    field from the standard elliptic-integral circular-loop formulas,
    Gauss-Legendre stacked over the sheet;
  * coils == J = NI/A over the cross-section, same loop kernel;
  * force on the neighbor mover via the magnetic-charge model
    F = sum_faces sigma_m * B_ext integrated over the pole faces
    (sigma_m = +/- Br/mu0), which is exact for uniform M in an external field;
  * everything validated against axisymmetric FEMM runs of the isolated
    mover through contour-flux window averages at r = 15/19.05/23 mm
    (mo_lineintegral(0); point evaluation hangs under Wine).

The FEMM validation also supplies empirical correction factors the analytic
model cannot produce:
  * k_washer  (w_o3_m12_r): far-field effect of the 2 mm mover washer,
  * k_out     (st_w03): axisymmetric leakage reduction by the 0.3 mm shell,
  * k_sleeve  (opp24c3 + shield sleeve): sleeve leakage reduction, and the
    sleeve run sweeps mover z to measure the cogging a sleeve introduces.

Run inside `nix develop` from the repo root, ALWAYS under `timeout`:
    python sim/qualify_crosstalk.py selftest
    python sim/qualify_crosstalk.py femm src opp24c3        # ~1 min each
    python sim/qualify_crosstalk.py femm src w_o3_m12_r
    python sim/qualify_crosstalk.py femm src st_w03_bare
    python sim/qualify_crosstalk.py femm src st_w03
    python sim/qualify_crosstalk.py femm sleeve             # 7 solves
    python sim/qualify_crosstalk.py compare                 # FEMM vs analytic
    python sim/qualify_crosstalk.py analyze                 # crosstalk report
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.special import ellipe, ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from femm import MU0, PROJECT_ROOT, RESULTS_DIR, z_path  # noqa: E402

PITCH = 19.05e-3          # m, adjacent actuator axis spacing
BR = 1.43                 # T, N52
MUR = 1.05
M_MAG = BR / MU0          # magnetization (charge model), A/m
K_SHEET = BR / (MU0 * MUR)  # equivalent solenoid sheet current = Hc, A/m
F_SCALE = 0.8             # N, force spec scale
F_FEEL = 0.010            # N, feel threshold
STROKE = 2.5e-3           # m, mover z_center +/- range


# ----------------------------------------------------------------- geometry
@dataclass
class Mag:
    a: float       # radius, m
    L: float       # length, m
    z_off: float   # center offset from mover reference, m
    sign: int      # +1 magnetized +z, -1 magnetized -z


@dataclass
class Cl:
    name: str
    r_in: float
    r_out: float
    z_bot: float
    z_top: float
    ni: float      # base ampere-turns at drive scale s=1


@dataclass
class Design:
    name: str
    magnets: list[Mag]
    coils: list[Cl]
    s08: float             # drive scale for 0.8 N in the worse direction
    sensor_z: float        # stator-fixed Hall position on own axis, m
    mover_bot_off: float   # bottom pole face offset from mover ref, m
    washer: tuple | None = None   # (r_out, z_bot, z_top) mover steel, FEMM only
    shell: tuple | None = None    # (r_in, r_out, z_half) stator steel, FEMM only
    note: str = ""


mm = 1e-3

DESIGNS: dict[str, Design] = {}


def _add(d: Design):
    DESIGNS[d.name] = d


# gains from results/opp24c3f_acd4.csv: worst (F(+1)-F(-1))/2 = 5.97 N/unit
_add(Design(
    "opp24c3",
    magnets=[Mag(4 * mm, 12 * mm, -6.1 * mm, +1),
             Mag(4 * mm, 12 * mm, +6.1 * mm, -1)],
    coils=[Cl("lo", 4.5 * mm, 7.75 * mm, -16.5 * mm, -8.0 * mm, -300.0),
           Cl("mid", 4.5 * mm, 7.75 * mm, -6.5 * mm, 6.5 * mm, +600.0),
           Cl("hi", 4.5 * mm, 7.75 * mm, 8.0 * mm, 16.5 * mm, -300.0)],
    s08=0.8 / 5.97,
    sensor_z=-15.5 * mm,       # on axis, below lowest mover-bottom (-14.6)
    mover_bot_off=-12.1 * mm,
    note="opposed pair D8x12, 3 air coils, pkg 35"))

# gains from results/w_o3_m12_r_wcgF.csv: worst (F(+1)-F(-1))/2 = 2.84 N/unit
_add(Design(
    "w_o3_m12_r",
    magnets=[Mag(4 * mm, 12 * mm, -7.0 * mm, +1),
             Mag(4 * mm, 12 * mm, +7.0 * mm, -1)],
    coils=[Cl("mid", 4.75 * mm, 7.75 * mm, -5.0 * mm, 5.0 * mm, +300.0),
           Cl("hi", 4.75 * mm, 7.75 * mm, 5.5 * mm, 15.5 * mm, -115.0),
           Cl("lo", 4.75 * mm, 7.75 * mm, -15.5 * mm, -5.5 * mm, -115.0)],
    s08=0.8 / 2.84,
    sensor_z=-16.0 * mm,       # below lowest mover-bottom (-15.5)
    mover_bot_off=-13.0 * mm,
    washer=(4 * mm, -1 * mm, 1 * mm),
    note="opposed pair D8x12 + 2 mm center washer, pkg 33"))

# gains from results/st_w03_fine_stloptfin.csv: worst g1 = 2.29 N/unit
_add(Design(
    "st_w03",
    magnets=[Mag(4 * mm, 12 * mm, 0.0, +1)],
    coils=[Cl("lo", 4.75 * mm, 7.65 * mm, -12.5 * mm, -0.5 * mm, -300.0),
           Cl("hi", 4.75 * mm, 7.65 * mm, 0.5 * mm, 12.5 * mm, +300.0)],
    s08=0.8 / 2.29,
    sensor_z=-9.5 * mm,        # in bore, below lowest mover-bottom (-8.5)
    mover_bot_off=-6.0 * mm,
    shell=(7.7 * mm, 8.0 * mm, 13.0 * mm),
    note="single D8x12 + 0.3 mm shell, pkg 28"))


# ------------------------------------------------------------- field kernels
def loop_B(a, rho, z):
    """B (T) of a circular loop radius a at origin carrying 1 A.

    Returns (B_rho, B_z); rho, z broadcastable arrays (m). Standard
    elliptic-integral form; rho -> 0 handled by the on-axis limit.
    """
    rho = np.asarray(rho, dtype=float)
    z = np.asarray(z, dtype=float)
    a = np.asarray(a, dtype=float)
    a, rho, z = np.broadcast_arrays(a, rho, z)
    Brho = np.zeros_like(rho)
    Bz = np.zeros_like(rho)
    on_ax = rho < 1e-9
    if np.any(on_ax):
        za, aa = z[on_ax], a[on_ax]
        Bz[on_ax] = MU0 * aa * aa / (2.0 * (aa * aa + za * za) ** 1.5)
    off = ~on_ax
    if np.any(off):
        r, zz, aa = rho[off], z[off], a[off]
        d2 = (aa + r) ** 2 + zz ** 2
        m = 4 * aa * r / d2                    # = k^2
        K, E = ellipk(m), ellipe(m)
        den = (aa - r) ** 2 + zz ** 2
        c = MU0 / (2 * np.pi * np.sqrt(d2))
        Bz[off] = c * (K + (aa * aa - r * r - zz * zz) / den * E)
        Brho[off] = c * zz / r * (-K + (aa * aa + r * r + zz * zz) / den * E)
    return Brho, Bz


def _gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def magnet_B(mag: Mag, z0, rho, z, n=32):
    """Field of one mover magnet with mover reference at z0. K = Hc sheet."""
    zc = z0 + mag.z_off
    zs, ws = _gl(n, zc - mag.L / 2, zc + mag.L / 2)
    rho = np.asarray(rho, dtype=float)[..., None]
    z = np.asarray(z, dtype=float)[..., None]
    br, bz = loop_B(mag.a, rho, z - zs)
    K = mag.sign * K_SHEET
    return (br * ws).sum(-1) * K, (bz * ws).sum(-1) * K


def mover_B(d: Design, z0, rho, z):
    Br = Bz = 0.0
    for mg in d.magnets:
        br, bz = magnet_B(mg, z0, rho, z)
        Br, Bz = Br + br, Bz + bz
    return Br, Bz


def coil_B(c: Cl, s, rho, z, nr=6, nz=16):
    """Field of one coil at drive scale s (total A-turns = s*ni)."""
    rs, wr = _gl(nr, c.r_in, c.r_out)
    zs, wz = _gl(nz, c.z_bot, c.z_top)
    J = s * c.ni / ((c.r_out - c.r_in) * (c.z_top - c.z_bot))  # A/m^2
    rho = np.asarray(rho, dtype=float)[..., None, None]
    z = np.asarray(z, dtype=float)[..., None, None]
    br, bz = loop_B(rs[:, None], rho, z - zs[None, :])
    w2 = wr[:, None] * wz[None, :] * J
    return (br * w2).sum((-1, -2)), (bz * w2).sum((-1, -2))


def coils_B(d: Design, s, rho, z):
    Br = Bz = 0.0
    for c in d.coils:
        br, bz = coil_B(c, s, rho, z)
        Br, Bz = Br + br, Bz + bz
    return Br, Bz


def B_xyz(field_rz, x, y, z):
    """Cylindrical source field -> Cartesian components at (x, y, z)."""
    rho = np.hypot(x, y)
    Br, Bz = field_rz(rho, z)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx = np.where(rho > 0, x / np.where(rho > 0, rho, 1.0), 1.0)
        cy = np.where(rho > 0, y / np.where(rho > 0, rho, 1.0), 0.0)
    return Br * cx, Br * cy, Bz


# ------------------------------------------------------- force (charge model)
def _face_quad(a, nr=8, nphi=16):
    rs, wr = _gl(nr, 0.0, a)
    phis = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    wphi = 2 * np.pi / nphi
    R, P = np.meshgrid(rs, phis, indexing="ij")
    W = (wr[:, None] * R) * wphi          # area weights r dr dphi
    return (R * np.cos(P)).ravel(), (R * np.sin(P)).ravel(), W.ravel()


def force_on_neighbor(field_rz, d_nb: Design, z_nb0, D=PITCH):
    """(Fx, Fz) on the neighbor mover (axis at x=D) in the source field.

    Charge model over each neighbor magnet's pole faces, sigma = sign*M on
    the top face, -sign*M on the bottom (M = Br/mu0). Neighbor washer/shell
    induced-moment forces neglected (second order in the ~mT external field).
    """
    xq, yq, wq = _face_quad(d_nb.magnets[0].a)
    Fx = Fz = 0.0
    for mg in d_nb.magnets:
        for zf, sgn in ((z_nb0 + mg.z_off + mg.L / 2, +mg.sign),
                        (z_nb0 + mg.z_off - mg.L / 2, -mg.sign)):
            bx, by, bz = B_xyz(field_rz, D + xq, yq, np.full_like(xq, zf))
            sigma = sgn * M_MAG
            Fx += sigma * np.sum(bx * wq)
            Fz += sigma * np.sum(bz * wq)
    return Fx, Fz


# ------------------------------------------------------------ flux (coupling)
def disk_flux(field_rz, xc, zc, a, nr=12, nphi=16):
    """Flux of the source field through a horizontal disk radius a centered
    at (xc, 0, zc) (axis-parallel normal). xc=0 -> own coil; xc=PITCH ->
    neighbor coil."""
    if xc == 0.0:
        rs, wr = _gl(nr, 0.0, a)
        _, bz = field_rz(rs, np.full_like(rs, zc))
        return float(np.sum(bz * rs * wr) * 2 * np.pi)
    xq, yq, wq = _face_quad(a, nr, nphi)
    _, _, bz = B_xyz(field_rz, xc + xq, yq, np.full_like(xq, zc))
    return float(np.sum(bz * wq))


def stack_linkage(d_src: Design, d_pick: Design, D):
    """Single-turn-pattern flux linkage of d_pick's coil stack in the field
    of d_src's coils driven at s=1. Turns weighted prop. to |NI| with winding
    sense sign(NI) (series-consistent). D=0 -> self inductance surrogate."""

    def f(rho, z):
        return coils_B(d_src, 1.0, rho, z)

    ni_max = max(abs(c.ni) for c in d_pick.coils)
    lam = 0.0
    for c in d_pick.coils:
        w = c.ni / ni_max
        zc = 0.5 * (c.z_bot + c.z_top)
        lam += w * disk_flux(f, D, zc, 0.5 * (c.r_in + c.r_out))
    return lam


# ------------------------------------------------------------------ FEMM I/O
HC_STR = f"{K_SHEET:.1f}"

# contour windows (mm): kind, then geometry
BR_RADII = (15.05, 19.05, 23.05)
BR_ZED = [(-24 + 4 * i, -20 + 4 * i) for i in range(12)]
BZ_ZC = (-16, -12, -8, -4, 0, 4, 8, 12, 16)
BZ_R = (18.05, 20.05)
BZA_ZC = (15, 18, 21, 24)
BZA_R = (0.05, 2.05)


def contour_defs(axis_windows=True):
    rows = []
    for rc in BR_RADII:
        for (z1, z2) in BR_ZED:
            rows.append(("Br", rc, z1, z2))
    for zc in BZ_ZC:
        rows.append(("Bz", zc, BZ_R[0], BZ_R[1]))
    if axis_windows:
        for zc in BZA_ZC:
            rows.append(("BzA", zc, BZA_R[0], BZA_R[1]))
    return rows


def _lua_header(L, out_csv):
    w = L.append
    w(f'LOG = "{z_path(out_csv)}.log"')
    w("function mark(s)")
    w('    local h = openfile(LOG, "a")')
    w('    write(h, s, "\\n")')
    w("    closefile(h)")
    w("end")
    w(f'handle = openfile("{z_path(out_csv)}", "w")')
    w('write(handle, "kind,c1,c2,c3,tot,avg\\n")')


class _Geo:
    """Segment-deduping geometry emitter (butted parts share edges)."""

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

    def poly(self, pts):
        n = len(pts)
        for i in range(n):
            (r1, z1), (r2, z2) = pts[i], pts[(i + 1) % n]
            self.seg(r1, z1, r2, z2)

    def label(self, r, z, mat, magdir=0, group=0, mesh=0.0):
        automesh = 1 if mesh <= 0 else 0
        self.w(f"mi_addblocklabel({r:.6g}, {z:.6g})")
        self.w(f"mi_selectlabel({r:.6g}, {z:.6g})")
        self.w(f'mi_setblockprop("{mat}", {automesh}, {mesh:.6g}, "", '
               f"{magdir}, {group}, 0)")
        self.w("mi_clearselected()")


def _emit_mover(g: _Geo, d: Design, z0_mm: float, with_shell: bool):
    for mg in d.magnets:
        zc = (z0_mm * mm + mg.z_off) / mm
        r_, l_ = mg.a / mm, mg.L / mm
        g.rect(0, zc - l_ / 2, r_, zc + l_ / 2)
        g.label(r_ / 2, zc, "Magnet", magdir=90 if mg.sign > 0 else 270,
                group=1, mesh=0.4)
    if d.washer is not None:
        r_o, zb, zt = (v / mm for v in d.washer)
        g.rect(0, z0_mm + zb, r_o, z0_mm + zt)
        g.label(r_o / 2, z0_mm + (zb + zt) / 2, "1018 Steel", group=1,
                mesh=0.3)
    if with_shell and d.shell is not None:
        r1, r2, zh = (v / mm for v in d.shell)
        g.rect(r1, -zh, r2, zh)
        g.label((r1 + r2) / 2, 0, "1018 Steel", group=0, mesh=0.25)


def _emit_contours(w, rows):
    for (kind, c1, c2, c3) in rows:
        if kind == "Br":
            p1, p2 = (c1, c2), (c1, c3)
        else:                                 # Bz / BzA: horizontal
            p1, p2 = (c2, c1), (c3, c1)
        w("mo_clearcontour()")
        w(f"mo_addcontour({p1[0]:.6g}, {p1[1]:.6g})")
        w(f"mo_addcontour({p2[0]:.6g}, {p2[1]:.6g})")
        w("tot, avg = mo_lineintegral(0)")
        w(f'write(handle, "{kind}", ",", {c1:.6g}, ",", {c2:.6g}, ",", '
          f'{c3:.6g}, ",", tot, ",", avg, "\\n")')
    w("mo_clearcontour()")


def gen_src_lua(d: Design, with_shell: bool, out_csv: str, fem: str) -> str:
    """Isolated mover (z0=0), far-field contour windows + mover Fz noise."""
    L: list[str] = []
    w = L.append
    _lua_header(L, out_csv)
    w("newdocument(0)")
    w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
    w('mi_getmaterial("Air")')
    w(f'mi_addmaterial("Magnet", {MUR}, {MUR}, {HC_STR}, '
      "0, 0.667, 0, 0, 1, 0, 0, 0)")
    if d.washer is not None or (with_shell and d.shell is not None):
        w('mi_getmaterial("1018 Steel")')
    g = _Geo(w)
    _emit_mover(g, d, 0.0, with_shell)
    # far-field mesh control: annulus covering the contour windows + a fine
    # axis box for the on-axis windows
    g.rect(14.5, -26, 24, 26)
    g.label(19.25, 25, "Air", mesh=1.0)
    g.rect(0, 14, 2.6, 25)
    g.label(1.3, 20, "Air", mesh=0.5)
    w("mi_makeABC(7, 60, 0, 0, 0)")
    g.label(30, 30, "Air")
    w('mark("built")')
    w(f'mi_saveas("{z_path(fem)}")')
    w('mark("analyze")')
    w("mi_analyze()")
    w('mark("solved")')
    w("mi_loadsolution()")
    w("mo_groupselectblock(1)")
    w("fz = mo_blockintegral(19)")
    w("mo_clearblock()")
    w('write(handle, "Fznoise,0,0,0,", fz, ",0\\n")')
    _emit_contours(w, contour_defs(axis_windows=True))
    w("mo_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


SLEEVE = (7.85, 8.15, 17.5)     # mm: r_in, r_out, z_half of candidate sleeve
SLEEVE_Z0 = [-2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5]


def gen_sleeve_lua(out_csv: str, fem: str) -> str:
    """opp24c3 mover + 0.3 mm 1018 sleeve (stator): cogging Fz(z0) with fine
    WST air, plus far-field contours at z0=0 for sleeve shielding."""
    d = DESIGNS["opp24c3"]
    L: list[str] = []
    w = L.append
    _lua_header(L, out_csv)
    for z0 in SLEEVE_Z0:
        w("newdocument(0)")
        w('mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)')
        w('mi_getmaterial("Air")')
        w(f'mi_addmaterial("Magnet", {MUR}, {MUR}, {HC_STR}, '
          "0, 0.667, 0, 0, 1, 0, 0, 0)")
        w('mi_getmaterial("1018 Steel")')
        g = _Geo(w)
        _emit_mover(g, d, z0, with_shell=False)
        r1, r2, zh = SLEEVE
        g.rect(r1, -zh, r2, zh)
        g.label((r1 + r2) / 2, 0, "1018 Steel", group=0, mesh=0.25)
        # fine WST air: U-shaped polygon around the mover (annulus + end
        # caps as ONE region; avoids coincident collinear segments)
        g.poly([(0, 14.7), (4.05, 14.7), (4.05, -14.7), (0, -14.7),
                (0, -17.4), (7.8, -17.4), (7.8, 17.4), (0, 17.4)])
        g.label(5.9, 0, "Air", mesh=0.4)
        g.rect(14.5, -26, 24, 26)
        g.label(19.25, 25, "Air", mesh=1.0)
        w("mi_makeABC(7, 60, 0, 0, 0)")
        g.label(30, 30, "Air")
        w(f'mark("z0={z0} built")')
        w(f'mi_saveas("{z_path(fem)}")')
        w(f'mark("z0={z0} analyze")')
        w("mi_analyze()")
        w('mark("solved")')
        w("mi_loadsolution()")
        w("mo_groupselectblock(1)")
        w("fz = mo_blockintegral(19)")
        w("mo_clearblock()")
        w(f'write(handle, "cog,{z0:.6g},0,0,", fz, ",0\\n")')
        if z0 == 0.0:
            _emit_contours(w, contour_defs(axis_windows=False))
        w("mo_close()")
        w("mi_close()")
    w("closefile(handle)")
    w("quit()")
    return "\n".join(L) + "\n"


def run_femm(lua_text: str, base: str, timeout: int = 1500):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    lua = base + ".lua"
    out_csv = base + ".csv"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua, "w") as f:
        f.write(lua_text)
    subprocess.run(["femm-lua", lua], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    return out_csv


FEMM_BASE = {
    "opp24c3": "xtalk_src_opp24c3",
    "w_o3_m12_r": "xtalk_src_w_o3_m12_r",
    "st_w03_bare": "xtalk_src_st_w03_bare",
    "st_w03": "xtalk_src_st_w03",
    "sleeve": "xtalk_sleeve_opp24c3",
}


def read_femm(tag: str):
    path = os.path.join(RESULTS_DIR, FEMM_BASE[tag] + ".csv")
    if not os.path.exists(path):
        return None
    out = []
    with open(path) as f:
        for row in csv.DictReader(f):
            out.append({"kind": row["kind"],
                        **{k: float(row[k]) for k in
                           ("c1", "c2", "c3", "tot", "avg")}})
    return out


# ------------------------------------------------- analytic window averages
def window_avg(field_rz, kind, c1, c2, c3, n=16):
    """Analytic surface-averaged normal flux density matching a FEMM contour
    window (all window coords in mm). Returns (flux_Wb, avg_T) with the
    normal = +r for Br windows and +z for Bz windows."""
    if kind == "Br":
        r = c1 * mm
        zs, ws = _gl(n, c2 * mm, c3 * mm)
        br, _ = field_rz(np.full_like(zs, r), zs)
        flux = 2 * np.pi * r * np.sum(br * ws)
        area = 2 * np.pi * r * (c3 - c2) * mm
    else:
        zc = c1 * mm
        rs, ws = _gl(n, c2 * mm, c3 * mm)
        _, bz = field_rz(rs, np.full_like(rs, zc))
        flux = 2 * np.pi * np.sum(bz * rs * ws)
        area = np.pi * ((c3 * mm) ** 2 - (c2 * mm) ** 2)
    return float(flux), float(flux / area)


def femm_vs_analytic(tag: str, field_rz, verbose=True):
    """Compare FEMM window averages (from tot/area, sign-robust) to analytic.
    Returns dict kind->(median |ratio|, list of (row, femm_avg, an_avg))."""
    rows = read_femm(tag)
    if rows is None:
        return None
    out = {}
    for r in rows:
        if r["kind"] in ("Fznoise", "cog"):
            continue
        if r["kind"] == "Br":
            area = 2 * np.pi * r["c1"] * mm * (r["c3"] - r["c2"]) * mm
            kind = "Br"
        else:
            area = np.pi * ((r["c3"] * mm) ** 2 - (r["c2"] * mm) ** 2)
            kind = r["kind"]
        femm_avg = r["tot"] / area
        _, an_avg = window_avg(field_rz, kind, r["c1"], r["c2"], r["c3"])
        out.setdefault(kind, []).append((r, femm_avg, an_avg))
    med = {}
    for kind, lst in out.items():
        ratios = [abs(f) / abs(a) for (_, f, a) in lst if abs(a) > 1e-8]
        med[kind] = float(np.median(ratios)) if ratios else float("nan")
        if verbose:
            print(f"  [{tag}] {kind}: median |femm/analytic| = {med[kind]:.3f}")
            for (r, f, a) in lst:
                print(f"    {kind} c=({r['c1']:g},{r['c2']:g},{r['c3']:g}) "
                      f"femm={f:+.4e} T  analytic={a:+.4e} T  "
                      f"ratio={f/a if abs(a) > 1e-9 else float('nan'):+.3f}")
    return med, out


def femm_ratio(tag_num: str, tag_den: str):
    """Window-by-window |avg| ratio of two FEMM runs (e.g. shell/bare)."""
    a, b = read_femm(tag_num), read_femm(tag_den)
    if a is None or b is None:
        return None
    key = lambda r: (r["kind"], r["c1"], r["c2"], r["c3"])
    bd = {key(r): r for r in b}
    ratios = {}
    for r in a:
        if r["kind"] in ("Fznoise", "cog"):
            continue
        rb = bd.get(key(r))
        if rb and abs(rb["tot"]) > 1e-12:
            ratios.setdefault(r["kind"], []).append(
                (r["c1"], r["c2"], r["c3"],
                 abs(r["tot"]) / abs(rb["tot"])))
    return ratios


# ------------------------------------------------------------------ analysis
def src_field(d: Design, z0):
    return lambda rho, z: mover_B(d, z0, rho, z)


def src_coil_field(d: Design, s):
    return lambda rho, z: coils_B(d, s, rho, z)


def scaled(field_rz, k):
    return lambda rho, z: tuple(k * b for b in field_rz(rho, z))


def analyze(k_src=None, s_in=(1.0, 3.0), sleeve_rows=None,
            k_sleeve=None, out_rows=None):
    """k_src: per-design source-field factor = median |FEMM/analytic| over
    the r=15-23 mm Br windows. Folds in (a) the ~+4% moment underestimate of
    the Hc-sheet model, (b) the washer (w_o3_m12_r), (c) the shell leakage
    reduction (st_w03, FEMM shell run vs bare analytic)."""
    k_src = k_src or {}
    D = PITCH
    zgrid = np.linspace(-STROKE, STROKE, 5)
    print(f"\n===== crosstalk analysis, pitch {D*1e3:.2f} mm =====")
    for name, d in DESIGNS.items():
        ksrc = k_src.get(name, 1.0)
        s_shield = s_in if name == "st_w03" else (1.0, 1.0)

        m_net = sum(mg.sign * BR * math.pi * mg.a ** 2 * mg.L / MU0
                    for mg in d.magnets)
        q_mom = sum(mg.sign * BR * math.pi * mg.a ** 2 * mg.L / MU0
                    * mg.z_off for mg in d.magnets)
        print(f"\n--- {name} ({d.note}) ---")
        print(f"  net dipole m = {m_net:+.3e} A m^2 ; "
              f"longitudinal quadrupole Q = sum(m_i z_i) = {q_mom:+.3e} A m^3")
        print(f"  source-side empirical factor k_src = {ksrc:.3f} ; "
              f"neighbor incoming shield S_in = {s_shield}")

        # ---------------- (a) field at neighbor sensor -------------------
        own = src_field(d, 0.0)
        _, bz_own = own(np.array([0.0]), np.array([d.sensor_z]))
        eps = 0.05 * mm
        _, bz_p = src_field(d, +eps)(np.array([0.0]), np.array([d.sensor_z]))
        _, bz_m = src_field(d, -eps)(np.array([0.0]), np.array([d.sensor_z]))
        slope = (bz_p[0] - bz_m[0]) / (2 * eps)      # T/m of mover motion
        gap = d.sensor_z - (d.mover_bot_off - STROKE)
        print(f"  sensor at z={d.sensor_z*1e3:+.1f} mm on own axis "
              f"(min gap to bottom pole {abs(gap)*1e3:.1f} mm):")
        print(f"    own field B_z = {bz_own[0]*1e3:+.1f} mT ; position gain "
              f"dB_z/dz_mover = {slope:.1f} T/m = {slope*1e3:.1f} mT/mm")

        bzx, bxx = [], []
        for z0 in zgrid:
            f = scaled(src_field(d, z0), ksrc)
            bx, by, bz = B_xyz(f, np.array([D]), np.array([0.0]),
                               np.array([d.sensor_z]))
            bzx.append(bz[0])
            bxx.append(bx[0])
        bzx, bxx = np.array(bzx), np.array(bxx)
        # neighbor sensor axial component ~unshielded even inside an open
        # tube; transverse component shielded by S_in
        dbz_stroke = bzx.max() - bzx.min()
        fmodc = scaled(src_coil_field(d, d.s08), ksrc)
        _, _, bz_c = B_xyz(fmodc, np.array([D]), np.array([0.0]),
                           np.array([d.sensor_z]))
        print(f"    crosstalk at neighbor sensor (axial B_z, unshielded by "
              f"neighbor tube): {bzx[len(bzx)//2]*1e6:+.1f} uT at src mid, "
              f"p-p over src stroke {dbz_stroke*1e6:.1f} uT")
        print(f"    crosstalk transverse B_x = {bxx[len(bxx)//2]*1e6:+.1f} uT "
              f"(/{s_shield[1]:.1f} if inside shielded bore)")
        print(f"    coil-drive crosstalk at sensor (src at 0.8 N drive): "
              f"B_z = {bz_c[0]*1e6:+.2f} uT (modulated at control bandwidth)")
        pos_err_um = abs(dbz_stroke / slope) * 1e6
        pos_err_coil_um = abs(2 * bz_c[0] / slope) * 1e6
        print(f"    -> position-sense corruption: {pos_err_um:.2f} um over "
              f"neighbor full stroke; {pos_err_coil_um:.3f} um from "
              f"neighbor 0.8 N coil drive")

        # ---------------- (b) force injected on neighbor mover -----------
        worst = (0.0, 0.0, 0.0, 0.0)   # |Fz|, Fx at that point, z_src, z_nb
        worst_fx = 0.0
        for z_src in zgrid:
            f = scaled(src_field(d, z_src), ksrc)
            for z_nb in zgrid:
                fx, fz = force_on_neighbor(f, d, z_nb, D)
                if abs(fz) > worst[0]:
                    worst = (abs(fz), fx, z_src, z_nb)
                worst_fx = max(worst_fx, abs(fx))
        div = s_shield
        print(f"  static magnet->magnet force on neighbor mover "
              f"(worst over +/-{STROKE*1e3:.1f} mm x2 grid):")
        print(f"    |F_z| = {worst[0]*1e3:.2f} mN (/{div[0]:.0f}..{div[1]:.0f} "
              f"incoming shield -> {worst[0]*1e3/div[1]:.2f}.."
              f"{worst[0]*1e3/div[0]:.2f} mN) at z_src={worst[2]*1e3:+.1f}, "
              f"z_nb={worst[3]*1e3:+.1f} mm")
        print(f"    |F_x| (lateral, rail load) up to {worst_fx*1e3:.1f} mN "
              f"(/shield same)")
        print(f"    vs 0.8 N drive scale: {worst[0]/F_SCALE*100:.2f}% ; "
              f"vs 10 mN feel threshold: {worst[0]/F_FEEL*100:.0f}%")

        # coil-drive force modulation on neighbor mover
        fmod_worst = 0.0
        fc = scaled(src_coil_field(d, d.s08), ksrc)
        for z_nb in zgrid:
            fx, fz = force_on_neighbor(fc, d, z_nb, D)
            fmod_worst = max(fmod_worst, abs(fz))
        print(f"    coil-drive modulation of neighbor force (src at 0.8 N): "
              f"|dF_z| = {fmod_worst*1e6:.1f} uN "
              f"({fmod_worst/F_FEEL*100:.3f}% of feel threshold)")

        # ---------------- coil-coil transformer coupling -----------------
        lam_self = stack_linkage(d, d, 0.0)
        lam_mut = stack_linkage(d, d, D)
        kk = abs(lam_mut / lam_self)
        print(f"  coil-coil coupling: lambda_self(1-turn pattern) = "
              f"{lam_self:.3e} Wb, lambda_mutual = {lam_mut:.3e} Wb "
              f"-> k = M/L ~ {kk:.2e}")

        # mover-motion induced EMF in neighbor coils (single-turn pattern)
        dz = 0.25 * mm
        def lam_from_mover(z0):
            f = scaled(src_field(d, z0), ksrc)
            ni_max = max(abs(c.ni) for c in d.coils)
            lam = 0.0
            for c in d.coils:
                zc = 0.5 * (c.z_bot + c.z_top)
                lam += (c.ni / ni_max) * disk_flux(
                    f, D, zc, 0.5 * (c.r_in + c.r_out))
            return lam
        dlam = (lam_from_mover(dz) - lam_from_mover(-dz)) / (2 * dz)
        v_key = 0.5   # m/s, fast keystroke
        print(f"  neighbor-keystroke EMF in coils: dlam/dz = {dlam:.2e} "
              f"Wb/m (1-turn) -> {abs(dlam)*v_key*1e6:.2f} uV/turn at "
              f"{v_key} m/s")

        if out_rows is not None:
            out_rows.append({
                "design": name, "k_src": ksrc,
                "m_net_Am2": m_net, "Q_Am3": q_mom,
                "B_own_sensor_mT": bz_own[0] * 1e3,
                "sense_gain_mT_per_mm": slope * 1e-3 * 1e3,
                "xtalk_Bz_sensor_uT": bzx[len(bzx)//2] * 1e6,
                "xtalk_Bz_pp_stroke_uT": dbz_stroke * 1e6,
                "xtalk_Bz_coil_uT": bz_c[0] * 1e6,
                "pos_err_stroke_um": pos_err_um,
                "pos_err_coil_um": pos_err_coil_um,
                "Fz_worst_mN": worst[0] * 1e3,
                "Fx_worst_mN": worst_fx * 1e3,
                "Fz_coilmod_uN": fmod_worst * 1e6,
                "k_coupling": kk,
            })

    # ---------------- (c) sleeve cogging -------------------------------
    if sleeve_rows:
        print("\n--- opp24c3 + 0.3 mm 1018 sleeve r7.85-8.15 z+/-17.5 "
              "(shield option) ---")
        cogs = {r["c1"]: r["tot"] for r in sleeve_rows if r["kind"] == "cog"}
        print("  cogging F_z(z0) [N] (raw, and antisym (F(z)-F(-z))/2):")
        for z0 in sorted(cogs):
            anti = ""
            if -z0 in cogs and z0 > 0:
                anti = f"   antisym {0.5*(cogs[z0]-cogs[-z0]):+.4f} N"
            print(f"    z0={z0:+4.1f} mm  Fz={cogs[z0]:+.4f} N{anti}")
        if k_sleeve:
            for kind, lst in k_sleeve.items():
                med = float(np.median([v[3] for v in lst]))
                print(f"  sleeve leakage ratio ({kind} windows, "
                      f"sleeve/bare): median {med:.3f}")


# ----------------------------------------------------------------- selftest
def selftest():
    # 1. thin solenoid center field vs smoke-test analytic (12.504 mT).
    # Dense z-quadrature: the loop kernel at the bore center is only ~a wide,
    # far narrower than in actual use (all real evaluations are >= 7 mm from
    # any winding).
    c = Cl("s", 4.75 * mm, 5.25 * mm, -50 * mm, 50 * mm, 1000.0)
    _, bz = coil_B(c, 1.0, np.array([0.0]), np.array([0.0]), nr=4, nz=256)
    ref = MU0 * (1000 / 0.1) * 0.05 / math.sqrt(0.05 ** 2 + 0.005 ** 2)
    print(f"solenoid Bz center: {bz[0]*1e3:.4f} mT vs {ref*1e3:.4f} mT "
          f"(err {abs(bz[0]/ref-1)*100:.2f}%)")
    assert abs(bz[0] / ref - 1) < 0.005

    # 1b. quadrature convergence at crosstalk distances: coil field at the
    # neighbor axis with default vs doubled node counts
    cc = DESIGNS["opp24c3"].coils[1]
    p_rho, p_z = np.array([PITCH]), np.array([-0.0155])
    b1 = coil_B(cc, 1.0, p_rho, p_z)[1][0]
    b2 = coil_B(cc, 1.0, p_rho, p_z, nr=12, nz=32)[1][0]
    print(f"coil quadrature at 19 mm: {b1*1e6:.4f} vs {b2*1e6:.4f} uT "
          f"(delta {abs(b1/b2-1)*100:.3f}%)")
    assert abs(b1 / b2 - 1) < 0.01
    mg0 = DESIGNS["opp24c3"].magnets[0]
    b1 = magnet_B(mg0, 0.0, p_rho, p_z)[1][0]
    b2 = magnet_B(mg0, 0.0, p_rho, p_z, n=96)[1][0]
    print(f"magnet quadrature at 19 mm: {b1*1e3:.5f} vs {b2*1e3:.5f} mT "
          f"(delta {abs(b1/b2-1)*100:.4f}%)")
    assert abs(b1 / b2 - 1) < 0.005

    # 2. magnet far field vs point dipole on axis at 100 mm
    mg = Mag(4 * mm, 12 * mm, 0.0, +1)
    _, bz = magnet_B(mg, 0.0, np.array([0.0]), np.array([0.1]))
    m_dip = BR * math.pi * mg.a ** 2 * mg.L / MU0
    ref = MU0 * 2 * m_dip / (4 * math.pi * 0.1 ** 3)
    # sheet model has m_eff = Hc*V = m/mu_r; compare against that
    ref_sheet = ref / MUR
    print(f"magnet on-axis 100 mm: {bz[0]*1e6:.3f} uT vs dipole(sheet) "
          f"{ref_sheet*1e6:.3f} uT (err {abs(bz[0]/ref_sheet-1)*100:.2f}%)")
    assert abs(bz[0] / ref_sheet - 1) < 0.02

    # 3. coaxial magnet-magnet force vs dipole-dipole at 100 mm (attractive
    # for like-oriented magnets -> F_z on the upper neighbor is negative)
    def f(rho, z):
        return magnet_B(mg, 0.0, rho, z)
    d_test = Design("t", [Mag(4 * mm, 12 * mm, 0.0, +1)], [], 0, 0, 0)
    fx, fz = force_on_neighbor(f, d_test, 0.100, D=0.0)
    m_src = m_dip / MUR          # sheet source
    m_nb = m_dip                 # charge-model neighbor uses Br/mu0
    ref = -3 * MU0 * m_src * m_nb / (2 * math.pi * 0.100 ** 4)
    print(f"coaxial force at 100 mm: {fz*1e3:.4f} mN vs dipole-dipole "
          f"{ref*1e3:.4f} mN (err {abs(fz/ref-1)*100:.2f}%), Fx={fx:.2e}")
    assert abs(fz / ref - 1) < 0.03
    print("selftest OK")


# ---------------------------------------------------------------------- cli
def main(argv):
    if not argv or argv[0] == "selftest":
        selftest()
        return
    if argv[0] == "femm":
        what = argv[1]
        if what == "src":
            name = argv[2]
            if name == "st_w03_bare":
                d, shell = DESIGNS["st_w03"], False
            else:
                d, shell = DESIGNS[name], True
            base = os.path.join(RESULTS_DIR, FEMM_BASE[name])
            run_femm(gen_src_lua(d, shell, base + ".csv", base + ".fem"),
                     base)
            print(f"done: {base}.csv")
        elif what == "sleeve":
            base = os.path.join(RESULTS_DIR, FEMM_BASE["sleeve"])
            run_femm(gen_sleeve_lua(base + ".csv", base + ".fem"), base)
            print(f"done: {base}.csv")
        return
    if argv[0] == "compare":
        print("== opp24c3 (bare opposed pair; validates analytic model) ==")
        femm_vs_analytic("opp24c3", src_field(DESIGNS["opp24c3"], 0.0))
        print("== st_w03_bare (single magnet; validates dipole case) ==")
        femm_vs_analytic("st_w03_bare", src_field(DESIGNS["st_w03"], 0.0))
        print("== w_o3_m12_r (washer effect = deviation from analytic) ==")
        femm_vs_analytic("w_o3_m12_r", src_field(DESIGNS["w_o3_m12_r"], 0.0))
        print("== st_w03 shell/bare leakage ratio ==")
        rr = femm_ratio("st_w03", "st_w03_bare")
        if rr:
            for kind, lst in rr.items():
                for (c1, c2, c3, v) in lst:
                    print(f"  {kind} ({c1:g},{c2:g},{c3:g}): "
                          f"shell/bare = {v:.3f}")
        return
    if argv[0] == "analyze":
        # pull empirical per-design factors from the FEMM validation runs
        k_src = {}
        for name, tag in (("opp24c3", "opp24c3"),
                          ("w_o3_m12_r", "w_o3_m12_r"),
                          ("st_w03", "st_w03")):
            r = femm_vs_analytic(tag, src_field(DESIGNS[name], 0.0),
                                 verbose=False)
            if r:
                k_src[name] = r[0].get("Br", 1.0)
        rr = femm_ratio("st_w03", "st_w03_bare")
        if rr and "Br" in rr:
            vals = [v[3] for v in rr["Br"]]
            print(f"st_w03 shell leakage factor (FEMM shell/bare, Br "
                  f"windows): median {float(np.median(vals)):.3f} "
                  f"(range {min(vals):.3f}..{max(vals):.3f})")
        print(f"per-design source factors k_src = "
              f"{ {k: round(v, 3) for k, v in k_src.items()} }")
        sleeve_rows = read_femm("sleeve")
        k_sleeve = femm_ratio("sleeve", "opp24c3")
        out_rows = []
        analyze(k_src=k_src, s_in=(1.0, 3.0),
                sleeve_rows=sleeve_rows, k_sleeve=k_sleeve,
                out_rows=out_rows)
        path = os.path.join(RESULTS_DIR, "crosstalk_summary.csv")
        with open(path, "w", newline="") as f:
            wcsv = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            wcsv.writeheader()
            wcsv.writerows(out_rows)
        print(f"\nwrote {path}")
        return
    raise SystemExit(f"unknown mode {argv}")


if __name__ == "__main__":
    main(sys.argv[1:])
