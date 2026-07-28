# Omniswitch actuator exploration — digest (2026-07-28)

Durable record of the FEMM exploration phase (axisymmetric magnetostatics,
`femm-lua` headless under Wine). Requirements are in `CLAUDE.md`; the three
finalists live as importable builders + metadata in **`sim/designs.py`**.

> **NOTE ON LOCATION**: `results/` is **gitignored** (generated FEMM outputs).
> This file is the durable record of the exploration — **move it to `docs/`**
> (e.g. `docs/EXPLORATION.md`) and commit it when you want it tracked.

All powers: ohmic W to hold 0.8 N in the *worse* direction anywhere in the
5 mm stroke, hot copper (2.1e-8 ohm-m), winding fill 0.6 unless noted.
"sym" = noise-cancelled symmetric power from the antisymmetric force gain;
"eval" = `analyze.evaluate` worst direction incl. cogging + F-I curvature.

## Finalists (all pass all specs on N52 at 20 C)

| # | design | topology | P_worst [W] (fill .6) | fill .5 | pkg [mm] | cogging [N] | tau | source |
|---|--------|----------|----------------------|---------|----------|-------------|-----|--------|
| 1 | `opp24c3` | opposed pair N52 D8x12 (N-N, 0.2 mm spacer), 3 air-core coils r 4.5-7.75, NI -300/+600/-300 | 0.36 (sym) - 0.41 (eval) | 0.44 | 35.0 (limit) | 0 (all-air) | 122 us | `sim/explore_aircore_dual.py`, `results/opp24c3f_acd4.csv` |
| 2 | `w_o3_m12_r` | opposed pair N52 D8x12 + 2 mm x 8 mm-OD 1018 washer between facing N poles (mover), ironless 3-coil r 4.75-7.75, NI +300/-115/-115 | 0.42 (sym) - 0.48 (eval) | 0.51 | 33.0 | 0 (no stator steel) | 92 us | `sim/explore_wildcard.py`, `results/w_o3_m12_r_wcgF.csv` |
| 3 | `st_w03` | single N52 D8x12, dual coil r 4.75-7.65, 0.3 mm 1018 shell r 7.7-8.0 = housing | 0.86 (sym) - 1.01 (eval) | 1.22 | 28.0 | 0.068 (real) | 162 us | `sim/explore_steel_shell.py`, `results/st_w03_fine_stloptfin.csv` |

Ranking logic: #1 lowest power + zero cogging but uses the whole 35 mm
budget; #2 nearly as efficient, 2 mm shorter, one extra part (washer);
#3 shortest/cheapest (shell is the housing), ~2.4x the hold power and a
small firmware-compensable cogging. Force-gain-vs-z tables for all three are
in `sim/designs.py` meta dicts (`gain_table`).

## Per-topology best numbers

Round-0 candidates (`sim/candidates.py`, tag `v0`):

| candidate | P_worst sym [W] | notes |
|---|---|---|
| `aircore_dual` (D8x12, 2 coils) | 1.07 | baseline |
| `aircore_single` (D8x12, 1 coil) | 4.43 | passes 5 W but marginal |
| `steel_shell_dual` (0.5 mm shell) | 0.83 | steel helps ~1.3x |
| `longmag_dual` (D8x20, pole coils) | 0.99 | long magnet helps |

Optimization rounds, best of each family:

| family (script) | best variant | P_worst [W] | pkg [mm] | cogging [N] | verdict |
|---|---|---|---|---|---|
| air-core dual, single magnet (`explore_aircore_dual`) | `a20` (D8x20, gap 4.5) | 0.70 sym | 34 | 0 | best non-opposed fallback |
| air-core opposed pair + 3 coils (`explore_aircore_dual`) | `opp24c3` | 0.36 sym | 35.0 | 0 | **FINALIST 1** |
| air-core single coil (`explore_aircore_single`) | `acs_best` = pc14m20c (Lc 14 @ zc 10.5, D8x20) | 1.61 sym | 32 | 0 | cheapest coil; keep as budget option |
| pole-washer wildcards (`explore_wildcard`) | `w_o3_m12_r` | 0.42 sym / 0.48 eval | 33.0 | 0 | **FINALIST 2** |
| plain pole washers on dual (`explore_wildcard`) | `w_washer_dual` | 0.95 sym | 27 | 0 | washers alone ~ +6% gain only |
| steel shell (`explore_steel_shell`) | `st_w03` (0.3 mm wall) | 0.86 sym / 1.01 eval | 28.0 | 0.068 | **FINALIST 3** |
| pot-core shell + end washers (`explore_steel_shell`) | `st_pot03/05/05L` | 0.74-0.80 sym | 29-33 | **0.40-0.88** | **REJECTED** (below) |

Other data points: `opp20c3` 0.48 W @ 31 mm (shorter opposed option),
`opp24w3` 0.40 W, `w_opposed3` (D8x10 + washer, equal NI) 0.57 W,
`w_o3_r` (ratio-tuned NI) 0.45 W, `st_w10` (1.0 mm wall) 0.88 W.

## Negative results (do not re-derive)

* **Pot-core end washers: REJECTED.** Adding steel end washers to the shell
  (st_pot03/05/05L) buys only 8-15 % power but adds a real, antisymmetric
  0.4-0.9 N cogging at the stroke ends (negative spring; st_pot03 ±0.88 N,
  st_pot05 ±0.77 N, st_pot05L ±0.40 N vs the 0.8 N full-scale force) —
  50-110 % of full scale is not firmware-compensable headroom.
* **FEMM harmonic solves are unreliable headless.** `explore_wildcard`'s AC
  eddy-check attempts hung at the first harmonic `mi_analyze()` in *every*
  configuration tried, including a coils+air-only problem (empty
  `results/w_o3_m12_r_wcgAC*.csv`, `results/acdbg_*`). `explore_steel_ac`'s
  linear-material harmonic sweep *did* run (`results/st_w03_stlac.csv`) —
  treat harmonic mode as fragile: keep materials linear, expect hangs, always
  use `timeout`, and keep the analytic diffusion bound as the fallback.
* **Eddy lag is a non-issue for the finalists.** st_w03 0.3 mm shell:
  harmonic flux-linkage check gives |dLam|/Lam < 0.03 % and phase < 0.4 deg
  at 500 Hz (equivalent first-order tau ~2 us); slitting the shell
  (zero-net-current parallel circuit) changes nothing measurable.
  w_o3_m12_r washer: analytic diffusion tau 0.15-0.3 ms at bias-realistic
  incremental mu_r 50-100 (3 ms even at unbiased mu_r 1000) and the primary
  force path (static PM field x coil current) bypasses the steel entirely —
  >30x margin on the 10 ms slew spec.
* **aircore_single family caps out ~1.6 W** (pole-centered pc14m20c) — fine
  vs the 5 W spec but 4x the opposed-pair power; keep only if a one-coil BOM
  is worth it.

## Verification caveats (bind all quoted numbers)

* **Stress-tensor noise floor.** Weighted-stress-tensor z-force on the mover
  carries a mesh-dependent constant offset: ~0.03-0.06 N typical with
  automesh (up to ~0.2 N on some solves), <= ~0.05 N with fine mesh
  (mesh_air 0.4 + 0.2-0.25 mm air in the mover-coil gap). The `ni_scale=0`
  rows report it directly. Any force claim below ~0.1 N (cogging!) needs
  fine mesh **and** antisymmetric +/-s differencing:
  gain = (F(+s) - F(-s))/2s cancels the offset exactly for air-core
  (F linear in I) and to first order for steel.
* **Fill factor.** `Coil.power()` uses fill 0.6, which is 10-20 % optimistic
  for layer-wound fine wire in a real bobbin — final powers are quoted at
  fill 0.5 too (x1.2). All finalists pass either way.
* **Boundary.** `abc_radius=60` converged: force shift <= 0.3 % vs R=120
  (checked for aircore_dual, opp24c3, w_o3_m12_r); R=40 shifts ~1.6 %.
* **Mesh.** aircore_dual convergence spot check (automesh vs 0.2 mm):
  force spread ~1.3 % at z=0, ~3.8 % at z=-2.5 across meshes; automesh vs
  finest within ~1 %. Percent-level force error is inherent — fine for
  ranking; the finalists' gain tables come from fine-mesh runs.
* **Analytic anchor.** `verify_dipole.py`: FEMM vs point-dipole force for a
  small magnet 10-20 mm from a current loop agrees to < 1 % after +/-NI
  differencing (raw ni=0 offset 3-7 mN even at 0.25 mm mesh).
* **Magnet grade/temperature.** All numbers are N52 at 20 C. Air-core power
  scales as 1/Br(T)^2 — presets and `power_scale()` in `sim/magnets.py`
  (e.g. N52 at 60 C costs ~1.13x power; N45SH ~1.24x vs N52 at 20 C).
* Winding turn count is a free parameter (R, L ~ N^2; P, tau invariant) —
  pick N last to match the driver and 24 V bus.

## File map

Simulation code (`sim/`):

| file | role |
|---|---|
| `femm.py` | core: Magnet/Coil/Steel/Model dataclasses, Lua codegen, headless runner, power model |
| `analyze.py` | spec compliance: worst-case power, cogging, tau, slew |
| `magnets.py` | magnet grade presets (Br(T), Hk(T), demag margin, power_scale) |
| `candidates.py` | round-0 candidate definitions (kept for history) |
| **`designs.py`** | **the three finalists: builders + verified metadata (start here)** |
| `explore_aircore_dual.py` | dual + opposed-pair air-core rounds (tags `acd1..4`); private stack/fine-sleeve Lua generator |
| `explore_aircore_single.py` | single-coil family, analytic prescreen + rounds (tags `acs_r1..3`, `acs_final`) |
| `explore_wildcard.py` | washer topologies; MoverModel multi-part mover (tags `wcg1..3`, `wcgF`); FEMM-harmonic-hang finding |
| `explore_steel_shell.py` | shell/pot-core family; ShellModel + fine annulus (tags `stlopt`, `stloptfin`) |
| `explore_steel_ac.py` | harmonic eddy-lag check via complex flux linkage (tag `stlac`) |
| `explore_meshconv.py`, `explore_abc_check.py`, `verify_dipole.py` | validation discipline: mesh, boundary, analytic anchor |

Key result CSVs (`results/`, gitignored; regenerate with the scripts above):

* Finalist reference sweeps: `opp24c3f_acd4.csv` (fine sleeve),
  `w_o3_m12_r_wcgF.csv` (fine, 11-pt z), `st_w03_fine_stloptfin.csv`
  (fine mesh + annulus, 11-pt z, 7 scales).
* Convergence/validation: `opp24c3r120_acd4.csv`, `opp24c3_acd4z.csv`,
  `w_o3_m12_r_wcgABC.csv`, `aircore_dual_meshconv_*.csv`,
  `aircore_dual_r{40,60,120}_abc.csv`, `verify_dipole_dipole.csv`,
  `st_w03_stlac.csv` (eddy), `acs_best_{fine,r120}_acs_final.csv`.
* Family sweeps: `*_v0.csv` (round 0), `a*_acd1.csv` / `*_acd2.csv` /
  `*_acd3.csv` (air-core dual rounds), `pc*/st*/mouth*_acs_r*.csv`
  (single-coil rounds), `w_*_wcg*.csv` (wildcards), `st_*_stlopt.csv`
  (shell family incl. the rejected pot cores).
* Empty-by-hang markers: `w_o3_m12_r_wcgAC.csv`, `w_o3_m12_r_wcgAC2.csv`,
  `acdbg_coils_wcdbg.csv` (harmonic mi_analyze hang evidence).

Fresh re-verification sweeps of the finalists write `results/<name>_designs.*`
(`python sim/designs.py <name> ...`, see its docstring).
