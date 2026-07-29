# omniswitch

A custom keyswitch that simulates **arbitrary force-distance curves** electronically:
a current-controlled solenoid (stator coil assembly) exerts force on a **permanent
magnet attached to the key**, with a displacement sensor closing the loop. Firmware
maps measured key position → commanded force in real time.

Current phase: **actuator design converged** — see `docs/DESIGN.md` (frozen
design: st_w10 — D8×12 N45SH magnet, dual coil, 1.0 mm steel shell/housing,
28×16 mm; builders in `sim/designs.py` + `sim/converge_shell.py`; exploration
digest in `docs/EXPLORATION.md`). Hall sensor placement frozen (cap-riding
in-bore stack — `docs/DESIGN.md` §Position sensor, study in
`sim/hall_position_study.py`). Next: Hall part selection, prototype, firmware.

## Hard requirements (actuator, per switch)

| Requirement | Value |
|---|---|
| Package OD | ≤ 16 mm (leaves inter-switch gap for cooling at 19.05 mm key pitch) |
| Full-force stroke | ≥ 0.8 N available everywhere across ≥ 5 mm of travel |
| Extra travel at lower force | bonus, not required |
| Package length incl. stroke | ≤ 35 mm |
| Continuous power at 0.8 N | ≤ 5 W at any point in the full-force stroke |
| Force slew | 0 → 0.8 N in < 10 ms from a bus voltage ≤ 24 V |
| Production cost | < $5/unit at ~100 units; easy assembly |
| Armature | **moving permanent magnet** (force ∝ current, bidirectional, zero cogging preferred) — NOT an iron plunger; precise open/closed-loop force control is the whole point |

Design implications worked out so far:
- Moving-magnet + coil gives F ∝ I (sign included). Electrical time constant
  τ = L/R is turn-count-invariant (~100 µs class for these sizes), so the 10 ms
  slew is easy **unless** solid steel is in the flux path — eddy currents in a
  solid return shell add a magnetic diffusion lag that can be ms-scale. Any
  steel-bearing variant must pass an AC (harmonic) force-vs-frequency check in
  FEMM. Thin-wall/slitted steel or ferrite are the fallbacks.
- Steel return shell ≈ doubles force per watt but adds a passive cogging force
  (magnet↔steel attraction at zero current) and mild F-vs-I nonlinearity. Cogging
  is firmware-compensable via the position sensor if ≪ 0.8 N — always report it
  (`ni_scale = 0` rows in sweeps).
- Winding turn count N only sets the voltage/current split: R ∝ N², L ∝ N²,
  required current ∝ 1/N, power and τ invariant. Pick N last, to match the
  motor driver and 24 V bus.

## Workflow

Commit AND push automatically at logical checkpoints (a completed study,
design decision, toolchain fix, doc update) — don't wait to be asked.

## Toolchain — everything through the Nix flake

All dependencies and builds are managed by `flake.nix` (this is a hard project
convention — do not install anything outside Nix).

- `nix develop` — dev shell: `femm-lua`, wine, xvfb-run, python3 (numpy/scipy/pandas/matplotlib).
- `nix build .#femm` — FEMM 4.2 Windows binaries, extracted reproducibly with
  innoextract (no Wine needed at build time). Unfree (Aladdin license), flake
  imports nixpkgs with `allowUnfree`.
- `femm-lua <script.lua>` — runs FEMM **headless** under Wine + xvfb.
  First run copies FEMM to `./.femm/app/` (FEMM writes state next to its exe)
  and creates a Wine prefix at `./.femm/wineprefix/` (~20 s, once). Both are
  gitignored.
- Flakes only see **git-tracked** files — `git add` new files before `nix build`.

## FEMM scripting conventions

- FEMM embeds **Lua 4.0** (not 5.x): stdlib is flattened into globals
  (`openfile`/`write`/`closefile`/`format`/`floor`…), no `%` operator, no
  `local function` niceties. **Don't hand-write model Lua** — generate it from
  Python (`sim/femm.py`); Python is the source of truth.
- Wine maps `Z:` → `/`; FEMM accepts forward slashes, so absolute paths are
  `"Z:" + posix_path` (`femm.z_path`).
- **Headless failure mode**: any FEMM Lua error raises a modal dialog → the run
  hangs forever. Always wrap runs in `timeout`, and use the generated
  `results/<run>.csv.log` step markers to identify the failing call after a hang.
- **APIs verified the hard way**: `mi_setcurrent` does NOT exist in FEMM 4.2
  (3.x leftover; calling it = modal dialog hang) — use
  `mi_modifycircprop("name", 1, amps)`. Point evaluation (`mo_getb`,
  `mo_getpointvalues`) hangs under Wine — use block integrals, contour
  integrals (`mo_lineintegral`), and `mo_getcircuitproperties` instead.
- The Bash-tool cwd persists across calls; femm-lua resolves its `.femm` state
  dir from `$PWD`. Run everything from the repo root with absolute paths, or
  you'll silently create stray `.femm`/`results` trees (this cost an hour once).
- Model conventions (see `sim/femm.py`): axisymmetric, mm, z = actuation axis;
  group 1 = mover (magnet + anything attached to key), group 0 = stator;
  force = weighted stress tensor z-force `mo_blockintegral(19)` (newtons);
  open boundary via `mi_makeABC(7, R, 0, 0, 0)` with R enclosing all mover
  positions; coils modeled as bulk copper regions driven by a 1-turn circuit
  carrying total ampere-turns NI; sweeps end with `quit()` or the headless run
  never exits.
- Winding power is computed analytically, not in FEMM:
  `P = ρ_cu · NI² · 2π·r_mean / (A_cross · fill)` with fill ≈ 0.6 and hot
  copper ρ = 2.1e-8 Ω·m (`Coil.power()`).
- Single-turn flux linkage from `mo_getcircuitproperties` is in the sweep CSVs;
  real-coil L = N² · Δλ/ΔI.
- Material library names verified in matlib.dat: `"Air"`, `"Copper"`,
  `"1018 Steel"`, `"NdFeB 52 MGOe"` (we usually use an explicit
  `mi_addmaterial` N52: Br 1.43 T, µr 1.05). AWG magnet-wire entries exist
  (`"30 AWG"` etc.) but aren't needed with the bulk-coil approach.

## Layout

- `flake.nix` — all deps/builds; packages: `femm`, `femm-lua`; dev shell.
- `sim/femm.py` — parametric model dataclasses (`Magnet`, `Coil`, `Steel`,
  `Model`), Lua codegen, headless runner (`run_sweep`), power model.
- `sim/*.lua` — hand-written smoke tests only (`hello.lua`,
  `smoke_solenoid.lua` — validates FEMM vs analytic finite-solenoid field).
- `results/` — generated .lua/.fem/.csv outputs (gitignored).

## Validation discipline

Never trust a new model class without: (1) analytic cross-check (smoke test),
(2) mesh-convergence spot check, (3) boundary-radius doubling check. Force
numbers from automesh + stress tensor can be percent-level off; that's fine for
ranking, tighten before committing to a design.
