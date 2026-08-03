# omniswitch actuator — frozen design record (2026-07-28)

Selected design: **st_w10** — single moving magnet, dual coil, 1.0 mm steel
flux-return shell doubling as housing and magnetic shield.

All numbers below are FEMM-verified (axisymmetric magnetostatics + harmonic
eddy study; pipeline validated to <1% against analytic models — see
docs/EXPLORATION.md and the sweep CSVs referenced in sim/converge_shell.py).

## Geometry (mm, radial from axis)

| Part | r | z | Notes |
|---|---|---|---|
| Magnet | 0 – 4.0 | ±6.0 (12 long), travels ±2.5 | D8×12 cylinder, axially magnetized, **N45SH** (N38UH if any 120 °C exposure) |
| Guide gap + bobbin | 4.0 – 4.75 | | PTFE/POM guide surface |
| Coils ×2 | 4.75 – 6.95 | lo −12.5..−0.5, hi 0.5..12.5 | opposed polarity, series, one continuous wire |
| Air | 6.95 – 7.0 | | |
| Shell = housing | 7.0 – 8.0 | ±13 | 1018 steel tube, 16 mm OD |
| Hall sensor asm | 0 – 3.5 | pedestal −13..−11.4, PCB −11.4..−10.6, package −10.6..−9.5 | on bottom cap, in bore; see §Position sensor |

Package: **28 × 16 mm** incl. 5 mm stroke and end caps (specs: ≤35, ≤16).

## Performance vs requirements

| Requirement | Result | Status |
|---|---|---|
| OD ≤ 16 mm | 16.0 | PASS (at limit) |
| ≥0.8 N everywhere on 5 mm stroke, both directions | gain 2.54–2.82 N per unit s (N52-sim); F ∝ I linear to <0.1%, no saturation to 2× drive | PASS, ~3× headroom |
| Length ≤ 35 mm | 28.0 | PASS |
| ≤5 W continuous at 0.8 N worst point | 1.53 W (20 °C) – 1.95 W (100 °C magnet), fill 0.5, N45SH, incl. cogging | PASS (thermal reality caps ~2–2.5 W continuous; we sit under it) |
| 0→0.8 N < 10 ms at ≤24 V | τ=137 µs, rise ≈49 µs + 5.7 µs shell eddy lag (measured, harmonic solve) | PASS ×~180 |
| <$5 @ 100 units | est. $3.5–4.5 (quotes pending) | PASS (verify) |
| Moving PM, precise force control | F ∝ I exact; cogging 0.161 N clean centering spring (odd in z) → firmware table | PASS |

## Why this one (decision trail)

- Air-core opposed-pair designs (opp24c3, w_o3_m12_r) were 2–4× more
  power-efficient BUT: (a) **crosstalk-fatal at 19.05 mm pitch** — the like-pole
  radial flux band injects 1.1–1.2 N on the neighbor's mover (quadrupole
  cancellation is a far-field story; 19 mm is deep near-field); (b) demag-stressed
  at the facing poles (opp24c3 unusable above ~70 °C even on N38UH); (c) opp24c3
  needs 31 N of permanent retention against pair repulsion (the washer variant
  cuts that to 2.4 N but keeps the crosstalk problem).
- Pot-core stator end washers: REJECTED — 0.4–0.9 N parasitic spring.
- 0.3 mm shell (st_w03): wall saturates at 1.85 T → only ~2× shielding;
  keyboard-matrix analysis says non-viable even with feedforward (residual up
  to 17% of 0.8 N).
- **st_w10's 1.0 mm wall stays at 1.21 T (unsaturated, µ_diff ≈ 325)**:
  measured source-leakage attenuation ×18 (worst neighbor force 25 mN,
  upper bound with unshielded victim), incoming transverse shielding S_in ≈ 23.
  Cost: +26% hold power vs thin shell — absorbed by margin.
- Demag: thicker return path slightly *raises* the magnet operating point.
  N52 disqualified in every design (margin ≤0.90). st_w10 margins (1.5×
  corner-peaking factor, bus-stall reverse overdrive s=1.27):
  N45SH 1.62/1.40/1.17 @ 80/100/120 °C; N38UH 2.20/1.94/1.65.
  Full record: results/demag_margins.csv (regenerate:
  `nix develop -c python sim/qualify_demag.py analyze`).
  Driver-overdrive addendum (2026-08-02): re-qualified at s=1.85 (555 A·t/coil
  = a 1.5 A driver on a 370 t/coil winding; covers any driver ≤1.5 A at that
  N): margins nearly unchanged — N45SH 1.58/1.36/1.14, N38UH 2.14/1.88/1.61
  @ 80/100/120 °C; N45SH stays ≥1.3 at 100 °C, so a 1.5 A-class driver is
  demag-safe. F(s=1.85) = 3.8 N worst pos/100 °C to 4.5 N center/80 °C,
  gain linear vs the stroke sweep to <0.2% (no shell saturation). Run:
  `converge_shell.py demag st_w10 1.85 demag185` → results/st_w10_demag185_*.

## Winding + drive (preliminary — re-quote after bobbin design)

~480 turns/coil AWG 33 heavy-build (fill 0.5), two coils series-opposed, one
continuous wire. R ≈ 27 Ω hot; hold at worst case ≈ 270 mA / 7.3 V (duty 0.30
on 24 V); bus-stall ≈ 0.9 A → stall force ≈ 2.6 N. Driver: DRV8871 for
bring-up, AT8870-class clone (~$0.3) for production; current sense via 0.3 Ω
shunt + INA181, mid-ON-time sampling; PWM 25 kHz phase-staggered in groups.
Force calibration: K(z) per-unit factory map on N45SH; temperature-compensate
K by −0.12 %/K of magnet temp using winding resistance as the thermometer.

## Thermal budget (2026-08-02, analytic + adversarially cross-checked)

Still air, standalone switch: shell sheds ~1.2–1.9 W at 60 K rise
(h_conv ≈ 9–10 W/m²K on 18 cm²; radiation 2 W/m²K bare steel → 7 painted —
paint/oxidize the shell if air-cooled); pins+PCB add a parallel ~30 K/W path
→ ~2–3 W total, confirming the 2–2.5 W "thermal reality" note. In-matrix the
3 mm inter-switch gaps ≈ boundary-layer thickness, so per-switch h drops and
the case exterior becomes the system bottleneck.

Water-cooled shell (jacket ~0.2 K/W, 25 °C coolant): the **winding's own
radial conduction binds**, and the magnet (no self-heating, floats to within
~6 K of the coil bore face) is the limiting node, not wire insulation.
R_wind = 1.10/k_eff K/W (exact uniform-generation annulus; dry round-wire
k_eff ≈ 0.15–0.35 W/mK, varnish/VPI-impregnated ≈ 0.6–0.96 — literature
range, the dominant uncertainty). Continuous capability at magnet ≤ 100 °C
(N45SH): as-built dry (0.05 mm air gap) ~12 W → 2.0 N; gap-potted only
~16 W → 2.3 N; **fully impregnated winding ~30–45 W → 3.2–3.8 N**, i.e. the
1.5 A driver ceiling (41–45 W hot) becomes available *continuously* — thermal
and driver limits coincide. Magnet PWM eddy heating: mW-class, negligible.
POM bobbin sits at magnet temp (~100 °C at the dry limit) → prefer PTFE-lined
or impregnated build. Impregnation is the single highest-leverage step; the
coil-shell gap may be potted freely (both stator, no moving interface).

## Keyboard integration (from the matrix feasibility study)

- 19.05 mm pitch is viable. Feedforward compensation of the neighbor force map
  is **mandatory** (firmware knows all key positions + currents); with st_w10
  leakage the residual is well under the 2–5%-of-0.8 N curve-shaping bar.
- Use **row-stripe magnet polarity** (alternate rows flipped): −40% adversarial
  / −81% deterministic neighbor map vs uniform. Checkerboard is worse — do not use.
- Sensor: neighbor field at own Hall ≈ 0.18 mT (3 µm apparent error) — negligible.
- Power: 104-key all-hold 0.8 N ≈ 62 W worst (24 V / 5 A supply with firmware
  global current budget + e-fuse; realistic typing ~5 W). ≥2 mF bulk +
  4.7 µF/key; 27–30 V TVS for release regen bursts.

## Shell slit: evaluated and REJECTED (2026-07-28)

Question: axial slit(s) in the shell to block eddy currents. Answer: **keep the
plain closed tube** — the slit removes almost nothing here, and costs real
performance. Studies: `sim/slit_shielding.py` (planar 2D transverse-shielding),
`sim/`+`results/st_w10_slitac.csv`, `st_w10_slitside.csv`.

- **The slit doesn't kill the eddy screening in this design.** Induced E is
  azimuthal, but because the coils are driven series-opposed the AC flux is odd
  in z and the net-winding circumferential mode carries ~zero current
  (FEMM: |I_net| ≤ 16 mA while distributed eddies dissipate 6.8 W at solve
  drive). The actual screening is zero-net-winding local loops (azimuthal out
  at one z, back at another, closing axially) — these survive a slit by closing
  along the slit edges. FEMM proof: a shell constrained to zero net current is
  digit-identical to the unslit shell at 20 Hz–50 kHz; a common-mode control
  drive (which this design never produces) shows the slit working (−81% lag).
  Analytic plate-mode estimate: 1 slit keeps ~65% of the lag, 2 slits ~37%.
- Benefit if slit anyway: 5.7 µs → ~2–4 µs lag (spec margin already 1750×) and
  1.4 mW → 0.3–0.7 mW PWM ripple loss (vs ~1.4 W hold). Noise.
- **Single-slit cost**: 27–54 mN constant-direction lateral preload on the mover
  (Maxwell stress of the missing wedge; kerf 0.2/0.4 mm) → 4–8 mN
  velocity-sign-dependent rail friction, 40–80% of the 10 mN feel threshold and
  NOT firmware-compensable (unlike cogging). Plus worst-orientation transverse
  shielding halves (S_eff 23 → 11–15).
- **Two symmetric slits**: side-load cancels to a 15–17 mN tolerance residual
  (±0.05 mm kerf, ±2°), BUT both circumferential shunt paths break —
  worst-orientation S_eff collapses to 4.0–5.6 (neighbor force 4.5–6.2 mN,
  approaching feel threshold before feedforward), the slit geometry is
  gap-reluctance-limited (high µ no longer helps), and the housing becomes two
  loose half-shells (~170× torsional derating).
- A slit aligned with the transverse field costs zero shielding — but a keyboard
  victim has orthogonal neighbors on both axes, so the worst orientation binds;
  45°-to-grid slit azimuth recovers ~25–40%.
- If a slit is ever *forced* (manufacturing): ONE slit, ≤0.2–0.3 mm kerf,
  bridged ends (one-piece tube), laser/EDM-cut DOM tube (never roll-formed
  open-seam as the OD-setting housing), azimuth 45° off-grid or toward board
  edge. If eddy screening ever genuinely matters, the fix is a ferrite return,
  not a slit.

## Position sensor: cap-riding in-bore Hall (added 2026-07-28)

Placement study: `sim/hall_position_study.py` → `results/hallpos.csv`
(15 magnetostatic solves over mover ±2.5 × drive {0,±1} × 9 sensor heights,
plus bare-magnet analytic check +2.2–3.7%, ABC-doubling ≤0.75% and
mesh-halving ≤0.44% on the own-signal channel; neighbor windows ≤4% at the
frozen in-bore heights, 5–18% at the deepest board-level rows — still far
inside the rejection margin). Board-level sensing (element under the package, e.g. z = −18)
**evaluated and REJECTED as primary**: below the shell the sensor loses the
×23 incoming shielding — neighbor-magnet pickup is ~85× the shielded in-bore
case (12-neighbor adversarial sum ≈ 1.6 mm apparent before feedforward,
breaching the 2–5% force-fidelity bar), and the own-coil term grows to
1.9 mm apparent. In-bore, every term is comfortable.

**Frozen setup — the sensor rides the bottom cap, inside the bore** (mm):
cap inner face −13.0 → molded pedestal Ø5.5 to −11.4 → sensor PCB Ø7.0 × 0.8
to −10.6 → SOT-23-class package (1.1) to −9.5, i.e. **1.0 clearance to the
magnet face at full press (−8.5)**. Sensing element sits 0.3–0.5 below the
package top → z ≈ −9.8…−10.0, matching the crosstalk study's assumed
z = −9.5 within 0.5 mm (its 0.18 mT / 3 µm neighbor numbers stand). The PCB
(Ø7.0) slides inside the Ø8.0 guide bore; 3–4 sensor pins + 2 coil tails exit
as through-cap pins (layout TBD) → the switch mounts through-hole and nothing
penetrates the steel shell.

At the element (N45SH scale): B ≈ 70–420 mT over the stroke, gradient
≈ 30–120 mT/mm → sub-2 µm noise floor at 50 µT-rms sensor noise. Own-coil
field ≈ 21 mT at rated drive (0.68 mm apparent) but is a clean k·I term
(constant over stroke to 1–2%, odd in drive sign to <0.2%) → calibrated
feedforward residual ~7–13 µm. Firmware obligations from the study: mid-PWM
synchronous Hall sampling (25 kHz ripple, ~0.3 mT class), rest-position
auto-zero (Hall offset drift 0.1–0.5 mT over temp), and the −0.12 %/K magnet
tempco must rescale the per-unit B(s) sensor map via the winding-resistance
thermometer (do NOT use sensor-IC built-in NdFeB compensation — the board/cap
sensor stays near ambient while the magnet heats). B(s) is a separate
factory-cal item from the K(z) force map (same rig).

**OPEN — range conflict at part selection**: no common linear Hall covers
420 mT (±300 mT is the ceiling: TMAG5170-A2 class). Resolve by (a) a
high-range part if one exists, or (b) shortening the pedestal — its height is
the single free parameter: h = 0 (PCB directly on the cap face) puts the
element at ≈ −11.5, B ≈ 40–237 mT (fits ±266/±300 with margin), clearance
2.6. Mechanical: the sensor stack must NEVER be the bottom travel stop —
stops (undesigned) must cap downward overtravel ≪ 1.0 beyond s = −2.5.
Fallback variant if through-hole mounting is dropped: flush PCB-mount with
the die in a cap cavity (≈ −13.4..−15 rows in hallpos.csv) — ~8× better
crosstalk than board-level, ~5× worse than in-bore.

## Open items before committing hardware

1. Vendor datasheet + quotes: D8×12 N45SH (margins near threshold use catalog
   Hk=0.9·Hci from memory — sim/magnets.py flags this), steel tube, winding.
2. Bobbin design → exact coil r_in; re-run winding numbers if r_in > 4.75.
3. Prototype validation: force-vs-current-vs-position rig; cogging table
   measurement; two-unit crosstalk measurement at 19.05 mm to confirm the ×18.
4. Hall part selection (placement now frozen — see §Position sensor). Must
   resolve the range conflict: element ≈ z −9.9 sees ~70–420 mT (N45SH),
   above the ±300 mT linear-Hall ceiling; either a high-range part or a
   shorter cap pedestal (h = 0 → element ≈ −11.5, 40–237 mT, clearance 2.6).
   Then: PWM-synchronous sampling, rest auto-zero, tempco-rescaled B(s) map.
5. End caps: steel caps would add axial shielding + a return path but also
   cogging — currently assumed polymer; simulate if sensor needs quieting.
   Bottom cap now carries the sensor pedestal + through-pin exits (§Position
   sensor); attachment to the shell and pin layout still undesigned.
