# CLAUDE.md — heliolinx-rhs developer guide

Guidance for Claude Code (and humans) working **inside this heliolinx checkout**.
This is the **`heliolinx-rhs`** fork (`github.com:rhstrauss/heliolinx-rhs`), a
performance/feature-extended version of the upstream `heliolinx/heliolinx`
HelioLinC3D C++ implementation. Active development happens on the **`omp_dev`**
branch.

For machine/pipeline/campaign context (Klone HPC, ATLAS runs, SLURM), see the
home-directory `~/CLAUDE.md` on each machine. This file is about the **code**.

## Build & install

```
cd src
make -j8            # build all PROGRAMS + libheliolinx.a
make install        # copies all PROGRAMS to ../bin  (PREFIX=.., override with PREFIX=/usr/local)
```

- Everything links against the single static lib `libheliolinx.a`
  (`solarsyst_dyn_geo01.cpp`, ~60k lines). Nearly every program `#include`s only
  `solarsyst_dyn_geo01.h` (the master header: all structs — `hldet`, `hlimage`,
  `tracklet`, `hlclust`, `shortclust`, `longpair`, `point3d`, configs — plus
  physical constants, time utils, Kepler/simplex tolerances).
- Generic Make pattern rules: adding a program = drop `src/foo.cpp` and add `foo`
  to the `PROGRAMS` line. Each program is `-std=c++11 -fopenmp -O3 -march=native`.
- Changing `solarsyst_dyn_geo01.h` rebuilds every program (all include it).

## Branch / multi-machine state

- Canonical branch: **`omp_dev`**. `main` tracks upstream; rebase/PR there only deliberately.
- Synced across **klone**, **gondor** (`rstrau@gondor.../heliolinx`, configured as a
  git remote on klone), and **origin** (GitHub `rhstrauss/heliolinx-rhs`). As of the
  last sync all three are at the same `omp_dev` tip. Sync path: klone ↔ origin ↔ gondor
  over git; prefer fast-forward / `-s ours` merges over force-push.
- **`bin/` is built from `src/` and is git-ignored.** Some SLURM pipelines invoke
  binaries straight out of `src/`, others out of `bin/`; after a `make` the `src/`
  binaries change immediately, so rebuilding mid-campaign can alter behavior of
  in-flight jobs. Be deliberate about when you `make`.

## Tool inventory (PROGRAMS)

`*_omp` = OpenMP-parallel variant. **(rhs)** = added/extended in this fork beyond
upstream.

### Tracklet creation
- **make_tracklets** — image-based pairing of detections within a night into
  pairs/tracklets; handles LSST/ATLAS column formats.
- **make_tracklets_omp** (rhs) — OpenMP-parallel make_tracklets (per nightly-chunk).
- **make_trailed_tracklets** — tracklets from *trailed* (streaked) detections of
  fast movers.
- **merge_det_catalogs** — combine multiple detection catalogs (multi-site) into one
  time-sorted 14-col catalog (parallel read + k-way heap merge).
- **merge_tracklet_files** — combine multiple tracklet file sets.
- **split_tracklets_by_time** — time-partition make_tracklets output (outim/pairdets/
  tracklets/trk2det) into ≤`-window`-day windows (`_split{NNN}`).
- **CROSS-OBSERVATORY linking (rhs, 2026-06):** `merge_pairs2`/`merge_trailpairs`
  take an optional `img_log`; when a tracklet candidate spans multiple observatories
  the GCR fit is redone with a per-image parallax correction over trial heliocentric
  distances, enabling cross-site links a single-site fit would reject. Helpers
  `gcr_fit_pairs`/`gcr_fit_trail`. Active in make_tracklets/make_trailed_tracklets.

### Linking — heliolinc proper
- **heliolinc** — main linker: iterate heliocentric (r, rdot) hypotheses, project
  tracklets to 3-D state vectors, k-d tree range-query clustering.
- **heliolinc_omp** — multithreaded heliolinc.
- **heliolinc_lowmem** — memory-efficient streaming variant (long windows / huge catalogs).
- **heliolinc_lowmem_omp** (rhs, primary workhorse) — OpenMP streaming: outer
  hypothesis loop parallel; each hyp writes its own `{outsum}_{N}.txt`/`{clust2det}_{N}.csv`.
  Flags `-streaming yes|no`, `-dedup yes|no`. **Cross-hypothesis dedup** is a
  parallel map-reduce funnel (`c0b6ad9`): per-thread `DedupSurvivor` maps + a
  deterministic tree-merge → byte-identical regardless of thread count, ~10× at 40
  threads. Keeps the **highest-metric** (lowest-RMS, best-fitting) representative of
  each identical detection set (`dedup_cluster_better`; flipped from lowest→highest
  2026-06 to match `link_dedup`/`link_dedup_funnel`).
- **heliolinc_estimate** (rhs) — predict RAM/runtime for a heliolinc_lowmem_omp run.
- **calc_heliohypmat** — build heliocentric hypothesis grids.
- **heliovane** — complementary algorithm for objects interior to Earth's orbit near
  90° phase.

### High-grading
- **helio_highgrade** — clustering-only pass over a hypothesis grid; emit every
  detection appearing in any cluster → smaller, cleaner catalog for a full run.
- **helio_highgrade_omp** (rhs) — parallel (atomic mark array). See `helio_highgrade_omp.md`.

### Post-processing (purify / planarity)
- **link_purify** — Keplerian orbit-fit every linkage; iteratively reject astrometric
  outliers to an RMS threshold; resolve overlapping linkages → non-overlapping final set.
- **link_purify_omp** (rhs) — parallel link_purify.
- **link_planarity** — link_purify + a fast coplanarity pre-screen (out-of-plane RMS)
  to reject outliers before the expensive orbit fit.
- **link_planarity_omp** (rhs) — parallel link_planarity; reads an `-lflist` of
  per-hyp pairs; merged output bit-identical to serial.
- **link_purify_chisq** (rhs) — link_purify variant using a chi-square cluster metric
  (`Hergetfit_vstar_chisq`). Includes the purify-glob **SIGSEGV fix** (guards
  `crossresid`/`alongresid` OOB after a failed orbit fit → mark as max-residual).
- **link_purify_chisq_omp** (rhs, 2026-06) — OpenMP link_purify_chisq. Phase A
  (per-cluster orbit fit) parallel; Phase B (greedy detection-conflict dedup) serial
  by nature. Phase-A exact-dup cull uses the parallel **`link_dedup_funnel`**.
  Optional **`-max_oop` planarity pre-cull**: when set (default off → identical to
  plain chisq), rejects non-coplanar linkages via `cluster_planarity_oop()` before
  fitting — fuses planarity + chisq in one pass. `-n_workers N`, `-heliovane`.
- **dedup_across_windows** (rhs) — cross-window idstring-based linkage deduplication.

### Output parsing / utilities
- **parse_clust2det** — expand cluster summaries into per-detection rows.
- **parse_clust2det_MPC80** — emit MPC 80-column astrometry.
- **modsplit_hlfile** — split/modify heliolinc output files.
- **parse_trk2det** — utility over tracklet-to-detection files.
- **label_hldet** — label detections in heliolinc output.
- **label_hldet_mpc** (rhs) — cross-match heliolinc detections against the MPC catalog.
- **analyze_linkage01a** — linkage diagnostics.

### Testing / injection
- **inject_fakes** (rhs, 2026-06) — synthetic small-body injection. Fully arg-driven
  (`-earth -obscode -incat -colformat -outim -orbits -population -outdets -truth ...`),
  no hardcoded paths → git-installable like any program. Generates fakes-only
  detection catalogs + a truth file, reusing the heliolinx ephemeris/projection. See
  `inject_fakes_design.md`.
- **tests/validate_pipeline.py** — end-to-end synthetic-MBA recovery test.

## Key library additions (in `solarsyst_dyn_geo01.cpp`)
- `link_dedup_funnel(inclust, inclust2det, out, out2det, nw)` — parallel tournament
  version of `link_dedup` (partition → leaf dedup → binary tree-merge; serial
  fallback nw≤1). Output identical to serial link_dedup (keeps highest-metric rep).
- `cluster_planarity_oop(...)` — thread-safe out-of-plane-RMS computation (extract of
  link_planarity geometry) for the `-max_oop` pre-cull.
- `gcr_fit_pairs` / `gcr_fit_trail` — factored GCR-fit helpers enabling the
  cross-observatory parallax tracklet feature.

## Recent major work (2026-05 → 2026-06, omp_dev)
1. `heliolinc_lowmem_omp` parallel cross-hypothesis funnel dedup (`c0b6ad9`).
2. Ported + improved `link_purify_chisq_omp` + `link_dedup_funnel` from gondor
   (`2fd48a1`); grafted the segfault fix; flipped heliolinc dedup to highest-metric.
3. `-max_oop` planarity pre-cull in `link_purify_chisq_omp` (`d157ead`).
4. `inject_fakes` made git-installable (`ec3d6ef`).
5. Cross-observatory parallax tracklet feature in make_tracklets (`1dc9a14`).
6. klone ↔ gondor ↔ GitHub reconciled to one `omp_dev` tip.

## Gotchas
- An experimental **H-magnitude filter** (`-minH/-maxH/-maxHspread` in
  `link_planarity_omp` + heliolinc config + `make_hypmask`/`prefilter_tracklets`) is
  **uncommitted WIP** in some working trees; not on `omp_dev` or synced.
- Pipelines may run binaries from `src/` (not `bin/`); a `make` changes those in place.
- The greedy detection-conflict dedup (link_planarity final pass, chisq Phase B) is
  inherently sequential — only the *exact-duplicate* dedup is funnel-parallelizable.
