# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

Heliolinx is a C++ implementation of the Heliolinc3D algorithm for asteroid discovery. It solves the "asteroid linking problem": identifying which telescope detections across multiple nights correspond to the same unknown asteroid. The core insight is projecting tracklets (intra-night detection pairs) to a common reference time under heliocentric distance hypotheses, then clustering the resulting 3D state vectors.

## Building

```bash
cd src
make -j4        # Build all executables
make install    # Install to ../bin (default)
PREFIX=/path make install  # Install to custom location
make clean      # Remove all binaries, objects, and .deps/
```

To build a single program (faster iteration during development):

```bash
cd src
make link_planarity_omp   # Only rebuilds this target and its dependencies
make install              # Copies all built programs to ../bin
```

Requires a C++11-capable compiler and OpenMP. The Makefile compiles with `-std=c++11 -fopenmp` and handles dependency tracking via `.deps/`. All programs link against `libheliolinx.a` (built from `solarsyst_dyn_geo01.cpp`), so changes to the shared library trigger recompilation of dependents.

## Building Python Bindings

```bash
pip install -e .       # Editable install
# or
python setup.py build_ext --inplace
```

The Python extension (`setup.py`) wraps `solarsyst_dyn_geo01.cpp` via PyBind11. The Python module source lives in `python/`. Platform-specific OpenMP flags are handled automatically (Darwin vs. Linux).

## Testing

Test data lives in a separate `heliolinx-aux/tests/` repository (not included here). The expected test inputs are:

- `test_TenObjects01a.csv` — detection catalog (10 simulated asteroids)
- `Earth1day2020s_02a.csv` — Earth ephemeris
- `ObsCodesNew.txt` — observatory codes
- `colformat_LSST_02.txt` — column format spec
- `heliohyp_rmb00a.txt` — heliocentric distance hypotheses

Canonical test sequence (from `README.md`):

```bash
make_tracklets -dets test_TenObjects01a.csv -earth Earth1day2020s_02a.csv \
  -obscode ObsCodesNew.txt -colformat colformat_LSST_02.txt

heliolinc -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
  -tracklets tracklets_TenObjects01a_01.csv -trk2det trk2det_TenObjects01a_01.csv \
  -obspos Earth1day2020s_02a.csv -heliodist heliohyp_rmb00a.txt -clustrad 2e5

link_purify -imgs outim_TenObjects01a_01.txt -pairdets pairdets_TenObjects01a_01.csv \
  -lflist TenObjects01a_lflist
```

Expected result: 10 linkages recovered (one per simulated asteroid).

## Architecture

### Pipeline Flow

```
Detection catalog
      ↓
make_tracklets  →  image list, paired detections, tracklets, trk2det
      ↓
[helio_highgrade]  →  pre-filtered tracklets (optional)
      ↓
heliolinc (or heliolinc_omp / heliolinc_lowmem / heliolinc_lowmem_omp)  →  cluster summary, clust2det
      ↓
link_purify / link_planarity  →  final deduplicated linkages
```

### Core Files

- **`src/solarsyst_dyn_geo01.h`** — The master header (~2100 lines). Defines all data structures (`hldet`, `hlimage`, `tracklet`, `hlclust`, `observatory`, `Point3D`, etc.), physical constants (AU, GM, speed of light), astronomical time utilities, and algorithm parameters (Kepler solver tolerances, simplex parameters, `NIGHTSTEP`, `INTEGERIZING_SCALEFAC`, etc.). Nearly every `.cpp` includes only this one header.

- **`src/solarsyst_dyn_geo01.cpp`** — Implements the bulk of the shared library (`libheliolinx.a`): coordinate transforms, orbit fitting (Herget method), k-d tree operations, ephemeris interpolation, and the heliocentric projection math (~56k lines).

- **`src/heliolinc.cpp`** — Main linking algorithm: iterates over heliocentric distance hypotheses, projects tracklets to 3D state vectors at reference time, runs k-d tree range queries to find clusters.

- **`src/make_tracklets.cpp`** — Pairs detections within a night into tracklets; handles LSST-style column formats.

- **`src/link_purify.cpp`** — Post-processing: orbit-fits each candidate linkage, rejects astrometric outliers, deduplicates overlapping clusters.

### Variants and Utilities

| Program | Purpose |
|---|---|
| `heliolinc` | Standard single-threaded version |
| `heliolinc_omp` | OpenMP parallelized version |
| `heliolinc_lowmem` | Reduced memory footprint |
| `heliolinc_lowmem_omp` | Low-memory + OpenMP parallelized version |
| `heliovane` | For interior asteroids at specific phase angles |
| `link_planarity` | Faster purification using coplanarity pre-screening |
| `helio_highgrade` | Pre-filters tracklets before full heliolinc run |
| `analyze_linkage01a` | Post-run linkage analysis |
| `calc_heliohypmat` | Calculates heliocentric hypothesis matrices |
| `label_hldet` | Labels heliolinc detections |
| `make_trailed_tracklets` | Makes tracklets from trailed sources |
| `merge_tracklet_files` | Merges multiple tracklet output files |
| `modsplit_hlfile` | Splits heliolinc output files |
| `parse_clust2det` | Parses cluster-to-detection mapping |
| `parse_clust2det_MPC80` | Parses clust2det in MPC 80-column format |
| `parse_trk2det` | Parses tracklet-to-detection mapping |
| `link_planarity_omp` | OpenMP parallelized planarity purification (production default) |
| `link_purify_chisq` | Chi-squared variant of link_purify |
| `label_hldet_mpc` | Cross-matches heliolinc detections against MPC catalog |
| `make_hypmask` | Creates hypothesis masks for targeted hypothesis runs |
| `prefilter_tracklets` | Pre-filters tracklets before heliolinc |
| `split_tracklets_by_time` | Splits tracklet files by time window |

### Key Concepts

- **Tracklet**: Two or more detections of the same object on the same night, defining a short arc.
- **Heliocentric hypothesis**: An assumed heliocentric distance (and its time derivative) used to project tracklets into 3D space.
- **State vector clustering**: After projection, nearby vectors in position+velocity space indicate the same real object.
- **`clustrad`**: The clustering radius parameter (in km) passed to `heliolinc`; controls sensitivity vs. false-positive rate.
- **`NIGHTSTEP`**: Minimum interval in days (defined in header) between detections to be counted as separate nights (currently 0.6 days).
- **`useunivar` modes**: Controls how heliolinc projects tracklets. Mode 0 = Keplerian, mode 2 = f-g function with `TrackletProjCache` (faster, uses heliocentric r, rdot, rdotdot triplets). Mode 2 is the production default.
- **`LinkPurifyConfig`**: Shared configuration struct (defined in the header) used by `link_planarity`, `link_planarity_omp`, `link_purify`, and `link_purify_chisq`. Holds all filtering thresholds (astrometric RMS, out-of-plane RMS, H-magnitude bounds, etc.).

## Critical: Internal Units

**Positions throughout the codebase are in km, not AU.** This includes `heliopos`, `image_log[].X/Y/Z` (observer positions), and all state vectors. The header defines `AU_KM = 1.495978700e8` for conversion. When computing physical quantities that need AU (e.g., absolute magnitude H = m - 5*log10(r_AU * delta_AU)), you **must** divide km values by `AU_KM`. Times are in MJD (days). Velocities are in km/sec or km/day depending on context — check carefully.

## Code Conventions

- All programs are self-contained `.cpp` files; each `main()` is in its own file.
- Shared functionality lives exclusively in `solarsyst_dyn_geo01.h/.cpp`.
- Programs parse their own command-line arguments by scanning `argv` for flag strings (e.g., `-dets`, `-earth`).
- Error handling follows a pattern of printing to `cerr` and returning nonzero, or (in newer code) skipping problematic clusters rather than exiting.
- Physical units throughout: **internal positions are always km** (see "Critical: Internal Units" above), times in MJD or Julian Days, angles in degrees or radians (check carefully). File I/O may use AU for hypothesis files (r, rdot) but internal computation is km.
