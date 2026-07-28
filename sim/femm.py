"""Parametric axisymmetric FEMM models for the omniswitch actuator.

Python is the source of truth: models are dataclasses, Lua scripts are
generated, run headless via the `femm-lua` wrapper (nix flake), and results
come back as CSV.

Conventions
-----------
* Units: mm for geometry, A-turns for coil excitation, N for force, W for power.
* Axisymmetric r-z plane; z is the actuation axis (key travel).
* Group 1 = mover (the magnet attached to the key); group 0 = stator.
* Force on the mover = weighted stress tensor z-force, mo_blockintegral(19).
* Coils are modeled as bulk-copper regions carrying total ampere-turns NI
  (circuit with 1 "turn"). Winding resistance/power is computed analytically
  in `Coil.power()` — it is independent of wire gauge for fixed NI and fill.
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
from dataclasses import dataclass, field

MU0 = 4e-7 * math.pi
RHO_CU_20C = 1.72e-8   # ohm*m
RHO_CU_HOT = 2.10e-8   # ohm*m, ~75 C winding — use for continuous-power checks

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")


def z_path(p: str) -> str:
    """POSIX path -> Wine path (Z: maps to /). FEMM accepts forward slashes."""
    return "Z:" + os.path.abspath(p)


@dataclass
class Coil:
    """Annular winding region, stator-fixed."""
    name: str
    r_in: float
    r_out: float
    z_bot: float
    z_top: float
    ni: float                # base ampere-turns (signed; sign sets polarity)
    fill: float = 0.6        # copper fill factor for layer-wound fine wire

    @property
    def area_mm2(self) -> float:
        return (self.r_out - self.r_in) * (self.z_top - self.z_bot)

    @property
    def r_mean_mm(self) -> float:
        return 0.5 * (self.r_in + self.r_out)

    def power(self, ni: float | None = None, rho: float = RHO_CU_HOT) -> float:
        """Ohmic power (W) to drive |NI| ampere-turns through this region.

        P = rho * NI^2 * 2*pi*r_mean / (A_cross * fill)
        """
        ni = self.ni if ni is None else ni
        a = self.area_mm2 * 1e-6
        rm = self.r_mean_mm * 1e-3
        return rho * ni * ni * 2 * math.pi * rm / (a * self.fill)


@dataclass
class Magnet:
    """Axially magnetized cylindrical magnet, the mover (group 1)."""
    radius: float
    length: float
    z_center: float          # position at sweep start
    br: float = 1.43         # N52
    mu_r: float = 1.05

    @property
    def hc(self) -> float:   # A/m
        return self.br / (self.mu_r * MU0)


@dataclass
class Steel:
    """Rectangular (in r-z) steel part. mover=True puts it in group 1."""
    name: str
    r_in: float
    r_out: float
    z_bot: float
    z_top: float
    mover: bool = False
    material: str = "1018 Steel"   # nonlinear BH from FEMM matlib


@dataclass
class Model:
    name: str
    magnet: Magnet
    coils: list[Coil]
    steels: list[Steel] = field(default_factory=list)
    abc_radius: float = 60.0       # open-boundary shell radius, must enclose all positions
    mesh_air: float = 0.0          # 0 = automesh; else max element size (mm) near mover
    fine_box: tuple | None = None  # (r_max, z_min, z_max) air box w/ mesh_air sizing

    def lua(self, z_positions: list[float], ni_scales: list[float],
            out_csv: str, fem_file: str) -> str:
        """Generate a Lua script sweeping mover position and excitation scale.

        The full geometry is regenerated per position (fresh document) instead
        of mi_movetranslate: FEMM re-meshes after a move anyway, and moving a
        group reliably requires assigning every node/segment to the group —
        an error-prone classic pitfall. Regeneration sidesteps it entirely.
        """
        m = self.magnet
        L: list[str] = []
        w = L.append

        # Step markers (append+close per mark so they survive a hang): FEMM
        # errors raise modal dialogs headless runs can't dismiss, so the .log
        # tail identifies the failing call.
        w(f'LOG = "{z_path(out_csv)}.log"')
        w("function mark(s)")
        w("    local h = openfile(LOG, \"a\")")
        w('    write(h, s, "\\n")')
        w("    closefile(h)")
        w("end")
        w(f'handle = openfile("{z_path(out_csv)}", "w")')
        flux_cols = "".join(f",flux_{c.name}" for c in self.coils)
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
            for s in self.steels:
                w(f'mi_getmaterial("{s.material}")')

            # Magnet (mover, group 1), magnetized +z (magdir 90 in the r-z plane)
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

            w(f"mi_makeABC(7, {self.abc_radius:.6g}, 0, 0, 0)")
            label(self.abc_radius * 0.5, self.abc_radius * 0.5, "Air")

            w(f'mark("z={z:.6g} built")')
            w(f'mi_saveas("{z_path(fem_file)}")')
            for s in ni_scales:
                for c in self.coils:
                    # NB: mi_setcurrent is a FEMM 3.x name and does NOT exist
                    # in 4.2 (calling it raises a modal dialog = headless hang)
                    w(f'mi_modifycircprop("{c.name}", 1, {c.ni * s:.6g})')
                w(f'mark("z={z:.6g} s={s:.6g} analyze")')
                w("mi_analyze()")
                w(f'mark("solved")')
                w("mi_loadsolution()")
                w("mo_groupselectblock(1)")
                w("fz = mo_blockintegral(19)")
                w(f'write(handle, {z:.6g}, ",", {s:.6g}, ",", fz)')
                for c in self.coils:
                    # single-turn flux linkage; real coil L = N^2 * dlambda/dI
                    w(f'ic, vc, lam = mo_getcircuitproperties("{c.name}")')
                    w('write(handle, ",", lam)')
                w('write(handle, "\\n")')
                w("mo_close()")
            w("mi_close()")
        w("closefile(handle)")
        w("quit()")
        return "\n".join(L) + "\n"


def run_sweep(model: Model, z_positions: list[float], ni_scales: list[float],
              tag: str, run_dir: str | None = None, timeout: int = 1800):
    """Generate, run headless, parse. Returns list of dicts (z, ni_scale, Fz)."""
    run_dir = run_dir or RESULTS_DIR
    os.makedirs(run_dir, exist_ok=True)
    base = os.path.join(run_dir, f"{model.name}_{tag}")
    lua_file, out_csv, fem_file = base + ".lua", base + ".csv", base + ".fem"
    if os.path.exists(out_csv):
        os.remove(out_csv)
    with open(lua_file, "w") as f:
        f.write(model.lua(z_positions, ni_scales, out_csv, fem_file))
    subprocess.run(["femm-lua", lua_file], check=True, timeout=timeout,
                   cwd=PROJECT_ROOT, capture_output=True)
    with open(out_csv) as f:
        return [{k: float(v) for k, v in row.items()}
                for row in csv.DictReader(f)]
