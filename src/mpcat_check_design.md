# mpcat_check / mpcat_index — design notes

`mpcat_check` answers one question fast: **which of these astrometric
measurements already exist in the Minor Planet Center observation archive?**
It is the heliolinx-native, git-installable C++ counterpart of the standalone
Python `mpcat_check` (`~/bin/mpcat_check`), and the rapid, pre-indexed analogue
of `label_hldet_mpc`.

## Why a new tool (vs. label_hldet_mpc)

`label_hldet_mpc` already does the right *match*: a 4-D (time + sky-position)
k-d tree cross-match of heliolinx detections against MPC 80-column observations,
gated by `-matchrad` (arcsec) AND `-timerad` (seconds), optionally requiring the
same observatory code. But it loads the obs80 files **fully into RAM**, which is
impossible for the ~40 GB / ~539 M-line MPC bulk archive (`NumObs.txt` +
`UnnObs.txt`).

`mpcat_check` keeps the identical match geometry but swaps the data path: it
queries a **pre-indexed, time-sorted binary catalog** (`mpcat.bin`, built once
by `mpcat_index`), mmap's it, and binary-searches **only the MJD window** spanned
by the input measurements. The k-d tree is built over that thin slice, so a query
touches a few million records, not 539 million.

## Two programs

### `mpcat_index` (build once)

```
mpcat_index -numobs NumObs.txt -unnobs UnnObs.txt [-extra CmtObs.txt ...] \
            -out mpcat.bin [-reserve 600000000]
```

- Streams each obs80 file line-by-line through `mpcat_parse_obs80line` (a quiet,
  fast parser added to `libheliolinx`), keeping only parseable optical RA/Dec
  observations. Two-line satellite/roving position lines, radar lines, and any
  malformed line are silently skipped (counted in the header).
- Holds all records in RAM, **sorts by MJD ascending**, and writes them as a raw
  `mpcdet[]` plus a text sidecar `mpcat.bin.hdr`.
- Memory-heavy: ~539 M records × `sizeof(mpcdet)`=56 B ≈ 30 GB resident. Run on a
  big-RAM node (see `KLONE_SCRIPTS/build_mpcat_index.sh`: `cpu-g2`, `--mem=64G`).
  `-reserve` pre-allocates the record vector to avoid a 2× reallocation spike.

### `mpcat_check` (per query)

```
mpcat_check -dets in.mpc80 -mpcat mpcat.bin \
            -matchrad 2.0 -timerad 5.0 -matchobscode 1 -out matched.csv \
            [-informat mpc80|pairdets] [-colformat colformat.txt] \
            [-unmatched_out new.csv]
```

- **Input** is either an MPC 80-column file or a heliolinx `pairdets` CSV.
  Auto-detected from the `-dets` extension (`.mpc80` → mpc80, else pairdets)
  unless `-informat` is given. The `pairdets` path reuses `read_detection_filemt2`
  and the same `-colformat` block as `label_hldet_mpc`.
- mmap's `mpcat.bin`, validates `sizeof(mpcdet)` against `mpcat.bin.hdr`
  (`recsize`), and `std::lower_bound`/`upper_bound`s the MJD-sorted array for
  `[min(MJD)-timerad, max(MJD)+timerad]` → a contiguous slice.
- Builds the 4-D k-d tree over the slice (same `DAY_TO_DEG_CONV=24` time scaling
  and `search_rad = sqrt(sky² + time²)` as `label_hldet_mpc`), then for each input
  detection range-queries and post-filters on `dt_sec ≤ timerad`,
  `sep_arcsec ≤ matchrad`, and (optionally) equal obscode, keeping the nearest
  qualifying MPC observation.

**Output** `-out` CSV, one row per input detection:
```
#MJD,RA,Dec,mag,obscode,idstring,in_mpc,mpc_packed,n_mpc_match,sep_arcsec,dt_sec
```
`in_mpc=1` rows carry the matched MPC packed designation, the number of matching
MPC observations, and the nearest match's separation (arcsec) and signed time
offset (sec); `in_mpc=0` rows are genuinely new. A stdout summary reports the
in-MPC / new split. `-unmatched_out` optionally writes just the new detections.

## `mpcat.bin` format

A flat array of fixed-layout `mpcdet` records (defined in
`solarsyst_dyn_geo01.h`), **sorted by MJD ascending**, with no header inside the
binary itself:

| field   | type            | source (obs80 cols) |
|---------|-----------------|---------------------|
| MJD     | double          | date, cols 16–32    |
| RA      | double (deg)    | cols 33–44 × 15     |
| Dec     | double (deg)    | cols 45–56 (signed) |
| mag     | float           | cols 66–70          |
| band    | char[5]         | col 71              |
| obscode | char[5]         | cols 78–80          |
| packed  | char[13]        | cols 1–12           |

`sizeof(mpcdet)` = 56 B (8-byte aligned). The sidecar `mpcat.bin.hdr` records
`nrec`, `recsize`, `tmin`, `tmax`, kept/skipped counts, and source files. The
record size is verified at query time; rebuild the catalog if `mpcdet` ever
changes layout.

## Catalogs on disk

- **klone:** `/gscratch/astro/rstrau/mpcat/{NumObs,UnnObs}.txt` → `mpcat.bin`.
- **gondor:** `/astro/store/shire/rstrau/mpcat/{NumObs,UnnObs}.txt` → `mpcat.bin`
  (also has `CmtObs.txt`; comets are out of scope for v1).

**Staleness:** the on-disk archives are snapshots (klone ~2024-03, gondor
~2025-11). Refresh with `download_mpcobs` and rebuild `mpcat.bin` when current
attribution matters. Do **not** rsync the ~30 GB `mpcat.bin` between machines —
regenerate it locally on each from that machine's obs80 files.

## Match tolerances

Defaults `-matchrad 2.0″`, `-timerad 5.0 s`, `-matchobscode 1` mirror
`label_hldet_mpc`. With the same observatory code required, a real re-detection
of an already-submitted observation matches to ≈0″/0 s (the input *is* that obs);
loosen `-matchrad`/drop `-matchobscode` only to find prior observations of the
same object at the same epoch from other sites.
