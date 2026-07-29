# reverse_design — design notes (heliolinc inversion / recoverability)

`reverse_design` answers the question heliolinc itself can only answer by
*running*: **given an asteroid with a known orbit, what heliocentric hypothesis
must be in the grid — and how finely must the grid be sampled — for heliolinc to
link it?** heliolinc is deterministic, so this is a closed-form inversion, not a
simulation. It is the *reverse* of `calc_heliohypmat` + the `trk2statevec`
projection: instead of enumerating a grid and discovering which objects fall out,
it takes one orbit and emits the exact `(HelioRad, R_dot, R_dubdot)` grid row that
recovers it, plus first-order **capture half-widths** giving the minimal grid
spacing that *guarantees* recovery.

It is the git-installable, heliolinx-native home for the analytic inversion first
developed and validated in the `NEO_discoverability` study (gondor); that study's
empirical validation and capture-scan campaign are summarized in §3–§4 below.

## What it is / how it installs

A standalone Python 3 CLI (no third-party deps; stdlib only). It is **not
compiled** — `make all` ignores it — but `make install` copies it to
`$(PREFIX)/bin/reverse_design` (the `.py` stripped, `chmod 755`), so it invokes
like any other heliolinx tool. Source of truth is `src/reverse_design.py`; `bin/`
is git-ignored. Physical constants (`K_GAUSS`, `AU_KM`, `SOLARDAY`) are kept
numerically identical to `solarsyst_dyn_geo01.h` so the inversion lands on the
same grid the C++ projection consumes.

## Usage

```
# Per-object prediction table (exact hypothesis + capture half-widths):
reverse_design -orbits orbits.txt -earth earth.csv -mjdref 60636.66819855 \
               -clustrad 2.0e5 -baseline_days 13.7858 -out predictions.csv

# Single-object validation probe: emit ONE grid row that should recover obj N:
reverse_design -orbits orbits.txt -earth earth.csv -mjdref 60636.66819855 \
               -emit_grid 0 obj0_hyp.txt
```

- `-orbits` : whitespace table `desig a e inc node argperi M epoch [H G]`
  (AU, deg, deg-mean-anomaly; `#` comment lines ok).
- `-earth`  : Horizons CSV. **Loaded only for interface parity** — the inversion
  is heliocentric and needs no Earth state; the flag is kept so the tool's
  signature matches the rest of the pipeline.
- `-mjdref` : the heliolinc reference epoch. For `-autorun 1` this is
  `0.5*MJDmin + 0.5*MJDmax` of the cadence (`outim_baseline()` helps compute it).
- `-clustrad` : heliolinc clustering radius, km (6-D L2). Default `2.0e5`.
- `-baseline_days` : full MJD span of the cadence; sets the half-baseline
  `T = baseline/2` used by the half-widths. Omitted → 1-day placeholder.

**Output** (`-out` CSV, one row per orbit):
```
#obj_index,desig,HelioRad_AU,R_dot_AUday,R_dubdot,
 dr_AU,drdot_AUday,dr_dubdot,
 min_grid_dr_AU,min_grid_drdot_AUday,min_grid_dr_dubdot
```
`HelioRad_AU,R_dot_AUday,R_dubdot` is the exact hypothesis; `dr_*,drdot_*,*dubdot`
are the per-axis capture half-widths; `min_grid_*` = `2×half-width` = the safe
maximum grid spacing per axis.

## 1. The hypothesis-grid column-3 convention (the crux)

A heliolinc hypothesis-grid row is `r(AU)  rdot(AU/day)  mean_accel`, mapped to
`hlradhyp{ HelioRad; R_dot; R_dubdot }` (`solarsyst_dyn_geo01.h` ~L636).

- **Col 1 `HelioRad`** = heliocentric distance r at MJDref, **AU**. Consumed as
  `heliodist = HelioRad * AU_KM` (`solarsyst_dyn_geo01.cpp` L39556).
- **Col 2 `R_dot`** = dr/dt at MJDref, **AU/day**. Consumed as
  `heliovel = R_dot * AU_KM` (km/day) — L39557.
- **Col 3 `mean_accel = R_dubdot`** = **normalized radial acceleration,
  dimensionless**, sign convention **"POSITIVE ACCELERATION IS INWARD TOWARD THE
  SUN."**

### Source proof — how it is WRITTEN (`calc_heliohypmat.cpp`)
```
L237:  g0 = GMSUN_KM3_SEC2/distkm/distkm;          // solar gravity, km/sec^2 (>0)
L252:  // WE TAKE POSITIVE ACCELERATION AS INWARD TOWARD THE SUN
L253:  maxacc = g0;   L254:  minacc = -g0;          // accel range is [-g0, +g0]
L277:  accelnorm = accelk/g0;                       // column 3 = accel / g0
L278:  outstream1 << dist << " " << velAU << " " << accelnorm << "\n";
```
For the rdot=0 reference rows the range is exactly `[-g0,+g0]` → column 3 ∈
[-1, +1], peaking at +1 for maximum inward acceleration. Matches `grid_master.txt`
(rdot≈0 rows sit near `mean_accel ≈ 1.0`).

### Source proof — how it is CONSUMED (`solarsyst_dyn_geo01.cpp` L39558)
```
helioacc = R_dubdot * ( -GMSUN_KM3_SEC2 * SOLARDAY^2 / heliodist^2 )
         = R_dubdot * ( -g0 )                       // g0 = GMsun/r^2 > 0
```
So the **physical outward radial acceleration** used in the projection is
`d²r/dt² = -R_dubdot · g0`. The explicit minus sign is the "positive = inward"
convention. heliolinc then projects each tracklet with the quadratic
`r_h(t) = r + rdot·dt + ½·helioacc·dt²` (`trk2statevec`, L33415).

### Inversion (orbit → exact hypothesis at MJDref)
1. Advance mean anomaly to MJDref; Keplerian elements → heliocentric ecliptic
   state (R, V) in AU, AU/day (`elements_to_rv`).
2. `HelioRad = |R|`.
3. `R_dot = (R·V)/|R|`.
4. `rddot = (|V|² − mu/r − R_dot²)/r`  (radial accel; from
   `r̈ = (|V|² + R·A − ṙ²)/r` with `A = −mu R/r³` ⇒ `R·A = −mu/r`).
5. `R_dubdot = −rddot / g0`,  `g0 = mu/r²`,  `mu = k²`, `k = 0.01720209895`.

A body on pure radial free-fall (`rddot = −g0`) ⇒ `R_dubdot = +1`, again
consistent with the reference rows. The ratio is dimensionless and
frame-independent, so it can be evaluated directly in AU units.

## 2. Capture half-widths (minimal guaranteed-recovery grid spacing)

heliolinc clusters 6-D state vectors — position (km) and velocity (km/s ×
`chartimescale`) — with a Euclidean L2 radius `clustrad` km at REF_GEODIST = 1 AU
(scaled ∝ geocentric distance). A hypothesis offset `(dr, drdot, dr_dubdot)`
perturbs the assumed heliocentric distance of a tracklet observed at time t by
```
δr_h(t) = dr + drdot·dt + ½·δ(helioacc)·dt² ,  δ(helioacc) = −g0·dr_dubdot ,  dt = t − tref.
```
Using the half-baseline `T = ½·(MJDmax − MJDmin)` as the characteristic |dt| and
requiring the edge tracklet to stay within `clustrad` of the reference epoch in
the position channel gives the first-order, single-channel half-widths:

| axis      | δr_h driver          | analytic half-width            |
|-----------|----------------------|--------------------------------|
| r         | `dr`                 | `clustrad_AU`                  |
| rdot      | `drdot·T`            | `clustrad_AU / T`              |
| r_dubdot  | `½·g0·dr_dubdot·T²`  | `2·clustrad_AU / (g0·T²)`      |

with `clustrad_AU = clustrad_km / AU_KM`; the implied minimal grid spacing is
`2 × half-width` per axis. For d050 (clustrad 2e5 km, T = 6.893 d):
`dr ≈ 1.34e-3 AU`, `drdot ≈ 1.94e-4 AU/day`, `dr_dubdot ≈ 0.198`. These are
deliberately **conservative** (a guaranteed-recovery lower bound): see §4.

## 3. Empirical validation — single-row probe

For each d050 object, emit ONLY its predicted 1-row hypothesis
(`reverse_design -emit_grid`), run `heliolinc_lowmem_omp -streaming no -dedup yes`
(clustrad 2e5, `-mingeodist 0.002 -useunivar 1 -minobsnights 2 -autorun 1`) →
`link_purify` → `parse_clust2det`, and check whether the object's `known_obj`
index appears in the parsed summary. **Result: 40 / 40 objects (100%) recovered by
their own single analytically-predicted hypothesis** — confirming both the
inversion math and the column-3 sign/normalization convention. Drivers + data live
in the `NEO_discoverability` study (`provenance/validate_reverse_design.sh`,
`data/inverse_d050.csv`).

## 4. Empirical capture widths vs. analytic (objects 0 and 7)

`capture_scan.sh` scans a 1-D offset ladder per axis and records the edge (the
smallest |offset| at which recovery fails).

| obj / axis     | analytic hw | empirical edge | ratio |
|----------------|-------------|----------------|-------|
| obj0 r         | 1.34e-3 AU  | 3.0e-2 AU      | ~22   |
| obj0 rdot      | 1.94e-4     | 8.0e-4         | ~4    |
| obj0 r_dubdot  | 0.198       | 0.60           | ~3    |
| obj7 r         | 1.34e-3 AU  | 1.2e-2 AU      | ~9    |
| obj7 rdot      | 1.94e-4     | 1.6e-3         | ~8    |
| obj7 r_dubdot  | 0.198       | 0.60           | ~3    |

The analytic half-widths are conservative by an axis-dependent factor; the
direction is correct and desirable (a guaranteed-recovery *lower bound*), so
`2×half-width` spacing safely guarantees recovery. Conservatism has two sources:
(1) `link_purify` re-fits a full Keplerian orbit and tolerates astrometric scatter
— recovery needs only `dbscan_npt` tracklets within `clustrad`, not all of them;
(2) the r-axis is most forgiving because a constant `dr` is *common-mode* (it
shifts every tracklet's projected distance nearly identically, absorbed by the
orbit re-fit, with extra slack from `clustrad ∝ geodist`); the true cluster-
breaking term for r is the *differential* shear, second order in `dt`, not the
first-order `dr` the table uses. rdot (linear in dt) and r_dubdot (quadratic in dt)
are genuinely differential, hence their much smaller factors (~3–8).

**Grid-design takeaway.** Use the analytic half-widths as a safe minimal spacing.
For a tighter (more efficient) grid, the r-axis spacing can be relaxed ~5–10× and
rdot/r_dubdot ~3–4× vs. the first-order estimate — but that trades the guarantee
for efficiency and must be re-checked per population/baseline.

## Scope & limits

- **Inversion ≠ detectability.** This predicts whether the *grid + clustrad* can
  cluster an object's tracklets; it does **not** model whether those tracklets
  exist (cadence, SNR, trailing, sky coverage). The observed ~10–15%
  "geometry-limited" recovery ceiling for the closest NEOs is that second factor,
  outside this tool.
- **Half-widths are single-channel, first-order.** They use only the *position*
  channel of the 6-D metric and treat the line-of-sight→heliocentric-distance map
  as ~1 (good for near-radial NEO geometry); the velocity channel and the
  geocentric-distance scaling of `clustrad` are not modeled — which is exactly why
  they are conservative. A tighter (velocity-channel + geo-scaled) model is the
  natural future extension.
- **Keplerian only.** The inversion uses two-body elements at MJDref; it matches
  heliolinc's own two-body projection but not perturbed dynamics over long arcs.

## Files

- `src/reverse_design.py` — the tool (inversion + half-widths + emit_grid).
- `NEO_discoverability/provenance/validate_reverse_design.sh` — 1-row recovery probe.
- `NEO_discoverability/provenance/capture_scan.sh` — per-axis capture-width scan.
- `NEO_discoverability/data/inverse_d050.csv`, `inverse_capture_d050.csv` — validation data.
