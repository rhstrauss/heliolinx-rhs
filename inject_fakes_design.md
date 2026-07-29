# `inject_fakes` — Synthetic Solar-System-Object Injection Tool: Scoping & Design

**Status:** Scoping (no code yet). Drafted 2026-05-22.
**Goal:** Measure heliolinx end-to-end recovery completeness for NEOs, MBAs (and, later, other SSOs) by injecting *physically motivated* synthetic objects into real ATLAS detection catalogs and scoring how many the pipeline recovers.

This plan was assembled from four parallel specialist investigations (small-body science/ADS literature, heliolinx format & CLI conventions, ATLAS data characterization, trailed-source modeling). The most important practical findings — and the conflicts between heliolinx's *designed* format and ATLAS's *actual* layout — are captured in §7.

---

## 1. Scope

### In scope (v1)
- **Catalog-level (detection-level) injection** — insert synthetic detection rows into the catalog that feeds `make_tracklets`. We are testing the *linking* pipeline, not pixel-level detection/photometry (ATLAS already solved that upstream).
- **Populations:** NEOs (Granvik et al. 2018 debiased model) and MBAs (Synthetic Solar System Model, Grav et al. 2011 / S3M).
- **Trailed and untrailed** synthetic detections.
- **Three output modes:**
  1. *Merged + labeled* — real detections + fakes, fakes tagged for recovery scoring.
  2. *Merged + unlabeled* — fakes indistinguishable from real detections (blind test); truth kept in a side-file only.
  3. *Fakes-only* — a catalog containing **only** injected objects, replaying the real catalog's cadence and per-exposure noise floor, with no stationary sources and no noise detections.
- **CLI that mirrors core heliolinx vocabulary** (`-indets`, `-colformat`, `-earth`, `-obscode`, `-outdets`, …).

### Deferred (post-v1)
- TNOs / Centaurs / Trojans (different heliocentric-distance + rate regime; needs a matching hypothesis grid). Clean extension hook only.
- Pixel-level injection (out of scope entirely — that tests a different pipeline).
- Split-detection modeling for ultra-fast movers (>100 deg/day trails crossing detection boundaries).

### Prior art on disk
`/gscratch/astro/rstrau/ATLAS/inject_test1/` (2026-03-28): `combined_inject.csv`, `hypotheses_inject.txt`, `truth_inject.csv`. A first manual attempt. Known defects to *not* repeat (see §7): `trail_len/trail_PA=0` (invalid pixel coords → trivially identifiable), `det_qual=1` (anomalous), hard magnitude cutoff instead of a detectability rollover. Also the in-repo `tests/validate_pipeline.py` (circular ecliptic orbits, no light-time) is a smoke test, **not** an adequate ephemeris engine for completeness work.

---

## 2. Architecture: a Python science core + reuse of `merge_det_catalogs`

Recommended split (resolves the "what language?" question — flag vocabulary mirrors heliolinx; the heavy science is where the libraries live):

```
            ┌─────────────────────────── inject_fakes (Python front-end) ───────────────────────────┐
            │                                                                                         │
  populations →  ① population sampler  →  ② ephemeris engine  →  ③ photometry+detectability  →  ④ catalog writer
  (Granvik,        (a,e,i,Ω,ω,M,H,G)       (pyoorb/OpenOrb;        (H-G mag, trailing loss,        (canonical 14-col
   S3M)                                     light-time,            per-exposure m5 cut,             + truth side-file)
                                            topocentric,           position scatter)
                                            per-obscode)                          │
                                                                                  ▼
                                          ⑤ MERGE: reuse existing `merge_det_catalogs` (C++)
                                             to interleave fakes + real into MJD-sorted output,
                                             format-preserving. Fakes-only mode skips the merge.
```

**Why reuse `merge_det_catalogs`:** it already does k-way MJD-sorted merging of multiple catalogs with per-catalog colformats into the canonical 14-column format, and is tested. The injection tool then becomes "generate a fakes catalog, hand it + the real catalog to merge_det_catalogs." Fakes-only mode is simply the Python-generated catalog with no merge step. This keeps format-preservation logic in one place and matches heliolinx idioms.

**Engine choices (from the science investigation):**
- **Ephemeris:** `pyoorb`/OpenOrb as production engine (returns topocentric RA/Dec, sky rates, phase angle, r, Δ, V for an MPC obscode in one vectorized call, **with iterated light-time**). `ASSIST`/`REBOUND` optional high-fidelity backend for long baselines. JPL Horizons (via `sbpy`/`astroquery`) as the validation oracle only (rate-limited; never per-object in bulk).
- **Photometry:** IAU H–G phase function (default G=0.15, optionally sampled); V→o/c color term matched to each exposure's band. `sbpy` for the phase-function math.
- **Populations:** Granvik 2018 residence-time maps (source-region weighting → (a,e,i) cell → H from source-specific slope; respect the q≲0.076 au super-catastrophic disruption depletion). MBAs drawn from S3M (gets Kirkwood gaps + family (a,e,i) structure for free). (Ω, ω, M) isotropic/uniform.

---

## 3. Module breakdown

| # | Module | Responsibility | Key inputs | Key outputs |
|---|--------|----------------|-----------|-------------|
| ① | `population` | Sample physically motivated orbits + H, G. NEO=Granvik 2018; MBA=S3M. | model files / S3M catalog; N, population, seed | table of (id, a,e,i,Ω,ω,M, epoch, H, G, taxon) |
| ② | `ephemeris` | For each orbit, compute topocentric RA/Dec, sky rate (μα·cosδ, μδ), α, r, Δ at each requested epoch/obscode. Iterated light-time. | orbits; exposure manifest (MJD, obscode, observer XYZ from outim / earth file) | per-(object,exposure) geometry rows |
| ③ | `photometry` | H–G apparent mag → band color term → trailing loss → per-exposure detectability rollover → drawn mag, mag_err, position scatter, trail shape. | geometry rows; per-exposure m5 + seeing + band | detected synthetic detections w/ all catalog fields |
| ④ | `catwriter` | Emit canonical 14-col catalog (+ matching colformat) and a truth side-file. Honor labeled/unlabeled. Enforce string-length limits. | detections | `fakes.csv`, `fakes_colformat.txt`, `truth.csv` |
| ⑤ | `merge` | (merged modes) call `merge_det_catalogs` on {real, fakes}. (fakes-only) no-op. | real catalog + fakes catalog | `outdets` |
| ⑥ | `score` | Post-pipeline: join recovered clusters back to truth → completeness vs H, a/e/i, rate, N-detections. | pipeline outputs + truth side-file | completeness tables/plots |

The **exposure manifest** (②) is the linchpin for both ephemeris epochs and fakes-only cadence replay. It comes straight from the `make_tracklets` **outim** file: one row per exposure = (MJD, boresight RA/Dec, obscode, observer XYZ/V, det index range, exptime=30s). Building `{(MJD,obscode)→outim_line}` lets us (a) place each fake in the right exposure group via the `image` field, and (b) decide which exposures actually saw each object.

---

## 4. CLI specification (mirrors heliolinx vocabulary)

```
inject_fakes
  -indets       <real detection catalog>        REQUIRED  (cf. make_tracklets -dets)
  -colformat    <colformat for -indets>         (omit if already canonical 14-col)
  -population   neo|mba|both                     REQUIRED
  -nobj         <N synthetic objects>            REQUIRED
  -popmodel     <Granvik/S3M model dir or file>  REQUIRED
  -earth        <heliocentric Earth ephemeris>   REQUIRED  (shared with make_tracklets)
  -obscode      <MPC observatory codes file>     REQUIRED  (shared with make_tracklets)
  -mjd          <reference MJD>                  (default: data midpoint, like -autorun)
  -outdets      <merged real+fake catalog>       OUTPUT (canonical 14-col)
  -outcolformat <colformat for -outdets>         OUTPUT (auto-written, always canonical)
  -outfakes     <fakes-only catalog>             OUTPUT (optional / fakes-only mode)
  -truth        <truth side-file>                OUTPUT (id↔orbit↔detection map)
  -labeled      yes|no   (default yes)           known_obj=fakeID vs known_obj=-1
  -trailed      yes|no|auto (default auto)        render trails when rate warrants
  -cadence_from_real yes|no (default yes)        replay outim cadence vs use orbit's own epochs
  -seed         <int>                            reproducibility
  -forcerun                                       push through non-fatal errors
  -verbose      <int>
```

Naming follows the suite idioms documented from the source: input data nouns (`-indets`/`-earth`/`-obscode`), `-out`-prefixed outputs, descriptor files (`-colformat`), and `yes|no` behavioral booleans in the style of `-streaming`/`-dedup`. Tool name follows `merge_det_catalogs` / `split_tracklets_by_time`.

---

## 5. Output modes in detail

**Canonical 14-column output** (always, regardless of input format), exact `merge_det_catalogs` layout:
```
#ID,MJD,RA,Dec,Mag,Band,ObsCode,trail_len,trail_PA,sigmag,sig_across,sig_along,known_obj,det_qual
fprintf: %s,%.7f,%.7f,%.7f,%.4f,%s,%s,%.2f,%.2f,%.4f,%.3f,%.3f,%ld,%ld
```

- **Merged + labeled** (default): `idstring=fake_<id>_<serial>`, `known_obj=<fakeID>`. Recovery scored by `pairdets.known_obj → clust2det → link_purify max_known_obj`.
- **Merged + unlabeled**: `known_obj=-1` (standard fill, indistinguishable). Identity recorded only in `-truth` side-file keyed by `idstring`. Blind efficiency test.
- **Fakes-only**: skip the merge. Output contains only injected detections, on exactly the real exposures (cadence + per-exposure m5 + seeing replayed from the real catalog), no stationary sources, no noise. Upper-bound recovery / ephemeris validation.

---

## 6. Recovery scoring (module ⑥)

Identity propagates: injected `known_obj`/`idstring` → preserved verbatim in `pairdets` (with `origindex`) → `clust2det` join → `link_purify` writes `max_known_obj` per cluster and the full per-detection block in `parse_clust2det` output. So:
- **Labeled:** a cluster is a recovery of object K if its detections carry `known_obj==K`; `max_known_obj` in the link_purify summary is a fast per-linkage flag.
- **Unlabeled:** join recovered cluster member `idstring`s against the truth side-file.

Report **completeness vs H, vs (a,e,i), vs apparent rate, and vs number of available detections** — the breakdowns where pipeline failure modes (tracklet formation, hypothesis-grid coverage, astrom-RMS cut) actually show up. Sanity floor: ATLAS-cadence completeness(H) should land near the published ATLAS debiased values (~88% at H<17.75) for an ATLAS-like configuration.

---

## 7. Critical correctness findings (read before implementing)

These are the non-obvious things the investigations surfaced. Several directly contradict the "clean" designed format.

1. **ATLAS `trail_len`/`trail_PA` are NOT trail geometry in the normal pipeline.** In the production (PSF) colformat, column 8 (`trail_len`) carries the **x pixel** coordinate and column 9 (`trail_PA`) the **y pixel** coordinate. Trail geometry (arcsec length, celestial PA) only appears in the *separate* trailed colformat (`ATLAS_colformat04trail.txt`). **The tool must be colformat-driven and populate whatever the target catalog's columns actually mean** — for normal-pipeline ATLAS that means plausible pixel coordinates (from the exposure WCS), not 0 and not arcsec. Setting them to 0 (as `inject_test1` did) makes fakes trivially identifiable.

2. **`sig_across`/`sig_along` do not affect linking.** `solarsyst_dyn_geo01.cpp` hardcodes `sigastrom=1.0` (with a "THIS IS CRUDE AND NEEDS FIXING" comment). The astrometric-uncertainty columns are cosmetic for recovery. **What matters is the actual position scatter we add to RA/Dec** — inject that directly and physically; the uncertainty columns just need to be plausible for output realism.

3. **Production `max_astrom_rms` is 0.5″, not the 2.0″ ATLAS default.** Current one-week scripts use 0.5″ (Pan-STARRS value). Real ATLAS linkages survive (median cluster RMS 0.41″), but the cut is tight — fakes with realistic position scatter near the limit may fail. **Test injections at both 0.5″ and 2.0″** and report which cut is binding.

4. **Use per-exposure limiting magnitude, not a constant.** Real ATLAS m5 swings 18.0–19.0 (o-band, T08) with weather/Vog/moon. Derive m5 per exposure from the real detections (95th-percentile mag of that image). This is the single biggest realism lever for completeness shape vs H.

5. **Detectability is a rollover, not a hard cut.** Use a logistic/erf in (m − m5) (width ~0.1–0.2 mag), or fit the Vereš et al. 2015 η_V curve to ATLAS. Optionally a rate-efficiency term η_U — but **do not double-count** it with trailing loss (one is a flux effect, the other a streak-detection effect).

6. **Trailing loss matters for the fast NEOs you care about.** `Δm = 1.25·log10(1 + a·x²/(1+b·x))`, `x = μ·t_exp/FWHM`. Detection (SNR) loss: a=0.67, b=1.16. Photometric loss: a=0.42, b=1.16 (Vereš et al. 2012, LSST-adopted). At ATLAS 30s, a 5°/day NEO trails ~6″ ≳ seeing.

7. **Trail astrometry is anisotropic.** Reported position = trail **midpoint** at exposure mid-MJD. Along-trail σ ≈ L/(2√3·SNR); cross-trail σ ≈ FWHM/(2√2·SNR). Inject scatter in the (along, cross) frame then rotate by PA to (RA, Dec) — not a single circular error. PA convention: E of N, mod 180° (trail axis has 180° ambiguity); verify "due-East motion → PA=90°".

8. **The `image` index is mandatory and must be correct.** It groups detections into the same exposure; a wrong/missing `image` prevents tracklet formation. Assign from the outim `{(MJD,obscode)→line}` lookup.

9. **String-length limits are asserted.** `idstring` ≤ 19 chars (`SHORTSTRINGLEN-1`), `band` ≤ 4, `obscode` ≤ 4 — violations abort in the hldet constructor.

10. **Ephemeris must be done properly.** Iterated light-time + topocentric (per-obscode parallax) are both arcsec-to-arcmin effects that will scramble or fabricate linkages. The current circular/no-light-time test is not usable here. Validate pyoorb vs Horizons to <0.1″ on a benchmark set before bulk runs.

---

## 8. Phased implementation plan

**Phase 0 — Environment & validation harness (small).** Stand up pyoorb/sbpy/astropy in the conda env. Validate pyoorb topocentric ephemerides vs JPL Horizons to <0.1″ for ~5 objects over a week. Confirm the exposure-manifest extraction from a real `outim` file.

**Phase 1 — Fakes-only generator (core science).** Modules ①②③④ for the *untrailed* case. Output a fakes-only catalog on a real one-week cadence. Validate: run it standalone through `make_tracklets → heliolinc_lowmem → link_planarity → link_purify` and confirm ~100% recovery of bright, well-sampled objects (this is the ephemeris/format correctness gate, replacing `validate_pipeline.py` with a physical engine).

**Phase 2 — Merge + labeling.** Wire module ⑤ (`merge_det_catalogs` reuse) and the labeled/unlabeled paths. Inject into a real highgraded one-week catalog; confirm format identity (make_tracklets runs unchanged) and that `known_obj` survives to `link_purify max_known_obj`.

**Phase 3 — Trailing model.** Add the trailed-source geometry/photometry/astrometry (§7.6–7.7) and the `auto` rate threshold. Decide trailed-colformat support (see §9). Calibrate trailing loss against a few known recovered ATLAS NEOs.

**Phase 4 — Recovery scoring + completeness product.** Module ⑥: completeness vs (H, a, e, i, rate, N-det). Run the full Granvik-NEO + S3M-MBA population against a representative window at both `max_astrom_rms` settings. Compare to the ATLAS debiased completeness as a sanity floor.

**Phase 5 — Hardening.** Reproducibility (seed), `-forcerun` semantics matching make_tracklets, docs (`inject_fakes.md`), and a regression test in `tests/`.

---

## 9. Open decisions for the user

1. **Implementation language for the science core** — recommend **Python** (pyoorb/sbpy/Granvik/S3M all live there) + **reuse C++ `merge_det_catalogs`** for the merge. Alternative: pure C++ (would mean re-implementing ephemerides — not recommended).
2. **v1 population scope** — recommend **NEO (Granvik 2018) + MBA (S3M)** only; defer TNO/Centaur. Confirm.
3. **Trailed colformat support in v1** — support both the normal (PSF, pixel-coord) and the dedicated trailed colformat, or normal-only first? Recommend normal-only in Phases 1–2, add trailed colformat in Phase 3.
4. **Which window(s)** to use as the injection testbed (e.g., `NEO_hg_onewk_01`, MJD 60630–60637, which has full 4-site coverage)?
5. **Where the tool lives** — `heliolinx/python/` + a thin CLI, or a standalone `src/` entry mirroring `merge_det_catalogs`?

---

## 10. References

- Granvik et al. 2018, Icarus 312, 181 (debiased NEO model) — arXiv:1804.10265
- Granvik et al. 2016, Nature 530, 303 (super-catastrophic disruption)
- Grav et al. 2011, PASP 123, 423 (S3M synthetic population)
- Denneau et al. 2013, PASP 125, 357 (Pan-STARRS MOPS injection/recovery) — arXiv:1302.7281
- Vereš et al. 2012, PASP 125, 1031 (trail fitting / trailing loss) — arXiv:1209.6106
- Vereš et al. 2015 (multi-survey detection efficiency incl. ATLAS) — arXiv:1511.07659
- Tonry et al. 2018 (ATLAS system) — arXiv:1802.00879
- ATLAS debiased NEO population 2024 — arXiv:2409.10453
- Holman et al. 2023 (ASSIST integrator) — arXiv:2303.16246
- LSST Solar System Yield 2025 — AJ, doi:10.3847/1538-3881/add685
- sbpy — sbpy.readthedocs.io; pyoorb/OpenOrb — github.com/oorb/oorb

*Internal: heliolinx `merge_det_catalogs`, `make_tracklets`, `split_tracklets_by_time`, `tests/validate_pipeline.py`; ATLAS catalogs under `/gscratch/astro/rstrau/ATLAS/`; prior art `inject_test1/`.*
