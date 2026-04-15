# helio_highgrade_omp vs helio_highgrade

`helio_highgrade_omp` is a drop-in OpenMP-parallel variant of `helio_highgrade`.
Both programs "high-grade" a detection catalog by keeping only those detections
that participate in at least one plausible heliocentric-motion cluster across a
supplied grid of radial-motion hypotheses, producing a much smaller input set
for the full linker. The only functional change is that `_omp` evaluates
hypotheses in parallel.

## What the two programs do

Given `-dets`, `-imgs`, `-tracklets`, `-trk2det`, `-obspos`, and a radial-motion
hypothesis file (`-radhyp`), each program:

1. Builds `heliodist / heliovel / helioacc` arrays from the hypothesis file.
2. For every hypothesis, converts the full tracklet set into six-element state
   vectors at the common reference MJD via one of
   `trk2statevec_fgfunc / trk2statevec_univar / trk2statevec_fgfuncRR /
   trk2statevec_univarRR` (selected by `-use_univar`).
3. Runs a KD-tree DBSCAN-style search (`highgrade_kdpairs`) on those state
   vectors to collect all detection indices that land in any cluster meeting
   the `-clustrad`, `-dbscan_npt`, `-mintimespan`, and geocentric-distance
   cuts.
4. Writes every detection so marked — de-duplicated across hypotheses — to
   `-outdets`.

## What differs

The CLI, input files, and output file format are identical. The differences
are internal and live in exactly two places.

### 1. Wrapper `src/helio_highgrade_omp.cpp`

- Header comment changed to describe the OMP variant.
- Adds a `-minobsnights` (aliases `-minnights`, `-minobsnite`) option to match
  the newer wrapper from the experimental branch. The serial program does not
  currently expose this flag.
- The final dispatch call changes from

      heliolinc_highgrade2(image_log, detvec, tracklets, trk2det,
                           radhyp, earthpos, config, minobsnum, outdets);

  to

      heliolinc_highgrade2_omp(image_log, detvec, tracklets, trk2det,
                               radhyp, earthpos, config, minobsnum, outdets);

  Everything else (argument parsing, input reading, output writing) is
  byte-identical to the serial wrapper.

### 2. Library function `heliolinc_highgrade2_omp` in `solarsyst_dyn_geo01.cpp`

This is a new function added next to the existing `heliolinc_highgrade2`. The
serial implementation is not touched.

The structural differences versus the serial `heliolinc_highgrade2` are:

- **Outer loop parallelism.** The hypothesis loop

      for(accelct=0; accelct<accelnum; accelct++) { ... }

  is replaced by

      #pragma omp parallel for schedule(dynamic) \
          default(none) shared(...)
      for(long ac=0; ac<accelnum; ac++) { ... }

  `schedule(dynamic)` handles the heavy hypothesis-to-hypothesis runtime
  variance seen with varying state-vector counts.

- **Per-thread scratch.** Inside the loop body, `allstatevecs` and
  `linkdet_temp` become local variables `local_statevecs` and `local_linkdet`
  so each thread owns its own memory and there is no write contention on
  them. `trk2statevec_*` and `highgrade_kdpairs` only read the shared input
  arrays and write to these thread-local vectors, so the parallel region is
  race-free.

- **Shared result via mark array, not growing vector.** The serial version
  maintains a `linkdet_indices` vector that it repeatedly extends, sorts, and
  deduplicates across hypotheses — an O(H·N log N) operation that also grows
  with every hypothesis. The OMP version replaces that entirely with

      vector<char> detmark(detnum, 0);
      ...
      #pragma omp atomic write
      detmark[idx] = 1;

  Each thread independently stamps a `1` into the mark array for every
  detection index returned by `highgrade_kdpairs`. After the parallel region
  the main thread walks `detmark` once and collects marked detections into
  `outdet`. This is O(H·k + N) with no per-hypothesis dedup pass and no
  accumulating vectors. Because each index is a distinct byte, the `#pragma
  omp atomic write` has no real contention — it only guarantees the store is
  not torn.

- **Logging is serialized.** The per-hypothesis "Hypothesis N: K detections
  found..." line, and the warning / error branches from `trk2statevec_*` and
  `highgrade_kdpairs`, are wrapped in `#pragma omp critical(hho_log)` so
  concurrent threads do not interleave partial lines on stdout/stderr.

- **Peak memory.** Bounded by `num_threads × (one hypothesis of state
  vectors + its highgrade_kdpairs scratch) + O(Ndet)` for the mark array,
  versus the serial version's `1 × (one hypothesis of state vectors) +
  O(Nkept_cumulative)`. For high-yield windows the OMP variant uses more
  memory during the parallel region, but it does not blow up with hypothesis
  count the way the serial accumulator could on very long runs.

- **Hypothesis ordering.** The serial version processes hypotheses strictly
  in input order. The OMP version processes them in whichever order threads
  pick them up. Because the final output is a deduplicated set keyed by
  detection index, the ordering does not matter for correctness, but log
  lines will interleave across hypotheses.

## What is unchanged

- The `HeliolincConfig`, `hlradhyp`, `hldet`, `tracklet`, and `longpair`
  structs, and every helper used inside the loop (`trk2statevec_fgfunc`,
  `trk2statevec_univar`, `trk2statevec_fgfuncRR`, `trk2statevec_univarRR`,
  `highgrade_kdpairs`, `earthpos01`).
- The serial `heliolinc_highgrade` and `heliolinc_highgrade2` library
  functions and the serial `helio_highgrade` binary.
- Output file format — the `-outdets` CSV is produced from `outdet` in the
  same way by both wrappers.

## Thread-count control

`helio_highgrade_omp` honors `OMP_NUM_THREADS` in the environment. It does
not currently expose a `-max_threads` CLI flag; set threads via

    OMP_NUM_THREADS=32 helio_highgrade_omp ...

as is already the convention for `heliolinc_lowmem_omp` and
`link_planarity_omp` on this branch.
