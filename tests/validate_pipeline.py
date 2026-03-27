#!/usr/bin/env python3
"""
validate_pipeline.py — End-to-end validation test for heliolinc.

Generates synthetic detections for N_ASTEROIDS simulated main-belt asteroids
on circular orbits near opposition, then runs the full pipeline:

    make_tracklets → heliolinc → link_purify

and verifies that every simulated object is recovered as a distinct linkage.

The synthetic data is self-consistent Keplerian mechanics (no precession,
no light-time correction, ecliptic-plane orbits) and is designed to exercise
every major stage of the pipeline, including the performance optimizations.

Usage:
    python3 tests/validate_pipeline.py [--keep-tmpdir] [--verbose]
"""

import argparse
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
AU_KM     = 149_597_870.7        # km per AU
GM_SUN    = 132_712_440_041.279  # km^3 / s^2
SOLARDAY  = 86_400.0             # seconds per day
MJDOFF    = 2_400_000.5          # JD - MJD

# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------
MJD_START = 59_200.0             # reference epoch (≈ 2021 Jan 06 UT)

# Four well-separated nights spanning >8 days (satisfies mintimespan=1 and
# minobsnights=3 with room to spare)
OBS_NIGHT_OFFSETS = [0, 3, 6, 10]          # days from MJD_START
DET_PAIR_GAP      = 0.5 / 24.0            # 30 minutes between the two detections

OBS_CODE = "500"                           # MPC geocentre; zero parallax offset

# Five main-belt asteroids at evenly spread heliocentric distances and
# slightly different starting angles so they appear at different RAs.
# theta0_deg is the heliocentric angle at MJD_START (measured in the
# ecliptic plane from the +x axis).  Setting it near 180° places each
# asteroid at approximate opposition.
ASTEROIDS = [
    {"id": "ast001", "r_h": 1.8, "theta0_deg": 181.0},
    {"id": "ast002", "r_h": 2.0, "theta0_deg": 185.0},
    {"id": "ast003", "r_h": 2.3, "theta0_deg": 173.0},
    {"id": "ast004", "r_h": 2.6, "theta0_deg": 177.0},
    {"id": "ast005", "r_h": 2.9, "theta0_deg": 169.0},
]
N_ASTEROIDS = len(ASTEROIDS)

# Tiny astrometric noise (arcsec → degrees), purely for realism;
# small enough that the test always passes.
ASTROMETRIC_NOISE_DEG = 0.1 / 3600.0

# Heliocentric hypothesis grid: circular-orbit hypotheses covering the full
# r_h range of the test asteroids.
HYPO_RH_MIN   = 1.5     # AU
HYPO_RH_MAX   = 3.2     # AU
HYPO_RH_STEP  = 0.1     # AU   (17 hypotheses — fast test)
HYPO_RDOT     = 0.0     # AU/day  (exactly circular)
HYPO_RDDOT    = 1.0     # dimensionless; =1 means exactly circular (GM/r^2)

# heliolinc parameters
CLUSTRAD      = 5e5     # km — generous radius; real runs use ~1-2e5
DBSCAN_NPT    = 3       # need 3 tracklets from 3 distinct nights
MINOBSNIGHTS  = 3

# ---------------------------------------------------------------------------
# Orbit mechanics helpers
# ---------------------------------------------------------------------------

def n_mean_motion(r_h_au):
    """Mean motion of a circular orbit at r_h_au, in rad/day."""
    n_earth = 2 * math.pi / 365.25          # rad/day
    return n_earth / (r_h_au ** 1.5)


def earth_state(mjd):
    """
    Heliocentric state of Earth at MJD *mjd*, assuming a circular orbit
    at 1 AU in the ecliptic plane (sufficient for a validation test).
    Returns (x, y, z [km], vx, vy, vz [km/s]).
    """
    n_e   = 2 * math.pi / 365.25
    phi   = n_e * (mjd - MJD_START)
    v_e   = math.sqrt(GM_SUN / AU_KM)      # ≈ 29.784 km/s
    x  =  AU_KM * math.cos(phi)
    y  =  AU_KM * math.sin(phi)
    vx = -v_e   * math.sin(phi)
    vy =  v_e   * math.cos(phi)
    return (x, y, 0.0, vx, vy, 0.0)


def asteroid_pos(mjd, r_h, theta0_deg):
    """
    Heliocentric position of an asteroid on a circular ecliptic orbit at
    heliocentric distance r_h AU, with initial angle theta0_deg at MJD_START.
    Returns (x, y, z) in km.
    """
    n_h   = n_mean_motion(r_h)
    theta = math.radians(theta0_deg) + n_h * (mjd - MJD_START)
    x = r_h * AU_KM * math.cos(theta)
    y = r_h * AU_KM * math.sin(theta)
    return (x, y, 0.0)


def geocentric_radec(ast_xyz, earth_xyz):
    """
    Topocentric RA / Dec (degrees) of an asteroid given heliocentric Cartesian
    positions.  For geocentre observations the topocentric = geocentric.
    """
    dx = ast_xyz[0] - earth_xyz[0]
    dy = ast_xyz[1] - earth_xyz[1]
    dz = ast_xyz[2] - earth_xyz[2]
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)
    ra  = math.degrees(math.atan2(dy, dx)) % 360.0
    dec = math.degrees(math.asin(dz / dist))
    return ra, dec

# ---------------------------------------------------------------------------
# File generators
# ---------------------------------------------------------------------------

def write_colformat(path):
    """Column-format descriptor for the detection catalogue."""
    with open(path, "w") as f:
        f.write("IDCOL 1\n")
        f.write("MJDCOL 2\n")
        f.write("RACOL 3\n")
        f.write("DECCOL 4\n")
        f.write("MAGCOL 5\n")
        f.write("BANDCOL 6\n")
        f.write("OBSCODECOL 7\n")


def write_obscodes(path):
    """
    MPC-style fixed-width observatory codes file.
    Only the geocentre (code 500) is needed.

    Parser reads:
      substr(0,3)  → 3-char code
      substr(4,9)  → longitude   (degrees E, 0 for geocentre)
      substr(13,8) → rho*cos(φ′) (= 0 for geocentre)
      substr(21,9) → rho*sin(φ′) (= 0 for geocentre)
    """
    header = "Code  Longitude   cos       sin       Observatory\n"
    # Carefully padded to put each value in the exact substring slice the
    # parser expects.  Total line length must be ≥ 30 characters.
    #
    # pos: 0-2   = "500"
    #      3     = " "
    #      4-12  = " 0.000000"  (9 chars)   → stod = 0.0  (longitude)
    #      13-20 = " 0.00000"   (8 chars)   → stod = 0.0  (plxcos)
    #      21-29 = " 0.000000"  (9 chars)   → stod = 0.0  (plxsin)
    #      30+   = "  Geocentre"
    line = "500  0.000000 0.000000 0.000000  Geocentre\n"
    # Verify the line is long enough and check the substring positions.
    assert len(line) >= 30, "ObsCodes line too short"
    with open(path, "w") as f:
        f.write(header)
        f.write(line)


def write_earth_ephemeris(path, mjd_start, mjd_end, step=1.0):
    """
    JPL-Horizons-style CSV ephemeris for Earth's heliocentric state.
    Covers [mjd_start - 5, mjd_end + 5] at *step*-day intervals so that
    heliolinc can interpolate anywhere within the observation window.

    Format expected by read_horizons_csv():
      header block (ignored)
      $$SOE
      JD, CalendarDate, X(km), Y(km), Z(km), VX(km/s), VY(km/s), VZ(km/s),
      ...
      $$EOE
    """
    with open(path, "w") as f:
        # Minimal header (read_horizons_csv skips everything before $$SOE)
        f.write("JPL Horizons — synthetic Earth ephemeris for heliolinc test\n")
        f.write("Columns: JDTDB, CalDate, X(km), Y(km), Z(km), VX, VY, VZ\n")
        f.write("$$SOE\n")
        mjd = mjd_start - 5.0
        while mjd <= mjd_end + 5.0:
            jd = mjd + MJDOFF
            x, y, z, vx, vy, vz = earth_state(mjd)
            # Calendar date placeholder (skipped by parser)
            caldate = "A.D. 2021-Jan-01 00:00:00.0000"
            f.write(f" {jd:.9f}, {caldate},"
                    f" {x:.6f}, {y:.6f}, {z:.6f},"
                    f" {vx:.9f}, {vy:.9f}, {vz:.9f},\n")
            mjd += step
        f.write("$$EOE\n")


def write_hypotheses(path):
    """
    Heliocentric distance hypothesis grid.
    All hypotheses use r_dot = 0 (circular orbit) and r_ddot = 1.0
    (the code interprets this as exactly one gravitational r̈ = -GM/r²).
    """
    r_h = HYPO_RH_MIN
    with open(path, "w") as f:
        while r_h <= HYPO_RH_MAX + 1e-9:
            f.write(f"{r_h:.4f} {HYPO_RDOT:.4f} {HYPO_RDDOT:.4f}\n")
            r_h += HYPO_RH_STEP


def write_detections(path, rng):
    """
    Detection catalogue CSV for all simulated asteroids.
    Each asteroid is observed on every night in OBS_NIGHT_OFFSETS, with two
    detections separated by DET_PAIR_GAP days (forming a tracklet).

    Columns: ID, MJD, RA, Dec, Mag, Band, ObsCode
    """
    rows = []
    det_id = 0
    for ast in ASTEROIDS:
        for night_offset in OBS_NIGHT_OFFSETS:
            for pair_idx in range(2):
                mjd  = MJD_START + night_offset + pair_idx * DET_PAIR_GAP
                ax, ay, az    = asteroid_pos(mjd, ast["r_h"], ast["theta0_deg"])
                ex, ey, ez, *_ = earth_state(mjd)
                ra, dec = geocentric_radec((ax, ay, az), (ex, ey, ez))
                # Add tiny Gaussian noise so detections aren't pixel-perfect
                ra  += rng.gauss(0, ASTROMETRIC_NOISE_DEG)
                dec += rng.gauss(0, ASTROMETRIC_NOISE_DEG)
                rows.append((ast["id"], mjd, ra, dec, 19.5, "r", OBS_CODE))
                det_id += 1

    # Sort by MJD then RA so make_tracklets sees a realistic time-ordered stream
    rows.sort(key=lambda r: (r[1], r[2]))

    with open(path, "w") as f:
        f.write("objID,MJD,RA,Dec,Mag,Band,ObsCode\n")
        for (oid, mjd, ra, dec, mag, band, obs) in rows:
            f.write(f"{oid},{mjd:.9f},{ra:.8f},{dec:.8f},{mag:.2f},{band},{obs}\n")

    return len(rows)


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def run(cmd, label, verbose=False, check=True):
    """Run a shell command, print output on failure (or always if verbose)."""
    if verbose:
        print(f"\n=== {label} ===")
        print("  " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or verbose:
        if result.stdout:
            for line in result.stdout.splitlines():
                print("  STDOUT:", line)
        if result.stderr:
            for line in result.stderr.splitlines():
                print("  STDERR:", line)
    if check and result.returncode != 0:
        print(f"\nFAIL: '{label}' exited with code {result.returncode}")
        sys.exit(1)
    return result


def count_data_lines(csv_path):
    """Count non-header, non-empty lines in a CSV."""
    n = 0
    with open(csv_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                n += 1
    # Subtract one for the header row if present
    with open(csv_path) as f:
        first = f.readline().strip()
    if first and not first.startswith("#") and not first[0].isdigit() and not first.startswith("-"):
        n -= 1
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--keep-tmpdir", action="store_true",
                        help="Do not delete the working directory after the test")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print full pipeline output")
    args = parser.parse_args()

    # Locate pipeline binaries relative to this script
    repo_root  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bin_dir    = os.path.join(repo_root, "src")
    make_trk   = os.path.join(bin_dir, "make_tracklets")
    heliolinc  = os.path.join(bin_dir, "heliolinc")
    link_pur   = os.path.join(bin_dir, "link_purify")

    for binary in (make_trk, heliolinc, link_pur):
        if not os.path.isfile(binary):
            print(f"ERROR: binary not found: {binary}")
            print("Please build first with: make -C src")
            sys.exit(2)

    # Deterministic RNG so the test is reproducible
    rng = random.Random(42)

    # Working directory
    tmpdir = tempfile.mkdtemp(prefix="heliolinc_test_")
    print(f"Working directory: {tmpdir}")
    if args.keep_tmpdir:
        print("(will be kept after test)")

    try:
        # ------------------------------------------------------------------ #
        # 1. Generate input files
        # ------------------------------------------------------------------ #
        det_file  = os.path.join(tmpdir, "detections.csv")
        earth_file = os.path.join(tmpdir, "earth_ephem.csv")
        obs_file   = os.path.join(tmpdir, "obscodes.txt")
        col_file   = os.path.join(tmpdir, "colformat.txt")
        hyp_file   = os.path.join(tmpdir, "hypotheses.txt")

        mjd_max = MJD_START + max(OBS_NIGHT_OFFSETS) + DET_PAIR_GAP
        write_colformat(col_file)
        write_obscodes(obs_file)
        write_earth_ephemeris(earth_file, MJD_START, mjd_max)
        write_hypotheses(hyp_file)
        n_dets = write_detections(det_file, rng)

        expected_tracklets = N_ASTEROIDS * len(OBS_NIGHT_OFFSETS)
        print(f"\nGenerated {n_dets} detections for {N_ASTEROIDS} asteroids "
              f"over {len(OBS_NIGHT_OFFSETS)} nights.")
        print(f"Expected tracklets: {expected_tracklets} "
              f"({N_ASTEROIDS} asteroids × {len(OBS_NIGHT_OFFSETS)} nights)")

        # ------------------------------------------------------------------ #
        # 2. make_tracklets
        # ------------------------------------------------------------------ #
        outim_file    = os.path.join(tmpdir, "images.txt")
        pairdets_file = os.path.join(tmpdir, "pairdets.csv")
        trk_file      = os.path.join(tmpdir, "tracklets.csv")
        trk2det_file  = os.path.join(tmpdir, "trk2det.csv")

        run([make_trk,
             "-dets",      det_file,
             "-earth",     earth_file,
             "-obscode",   obs_file,
             "-colformat", col_file,
             "-outim",     outim_file,
             "-pairdets",  pairdets_file,
             "-tracklets", trk_file,
             "-trk2det",   trk2det_file,
             "-forcerun",  "1"],
            label="make_tracklets", verbose=args.verbose)

        # Count tracklets (each line in the tracklets file is one tracklet)
        n_tracklets = count_data_lines(trk_file)
        print(f"\nmake_tracklets: found {n_tracklets} tracklets "
              f"(expected {expected_tracklets})")
        if n_tracklets < expected_tracklets:
            print(f"FAIL: expected at least {expected_tracklets} tracklets, "
                  f"got {n_tracklets}")
            sys.exit(1)

        # ------------------------------------------------------------------ #
        # 3. heliolinc
        # ------------------------------------------------------------------ #
        # Reference MJD near the midpoint of the observation window
        mjd_ref     = MJD_START + (max(OBS_NIGHT_OFFSETS) / 2.0)
        outsum_file  = os.path.join(tmpdir, "clusters.csv")
        c2d_file     = os.path.join(tmpdir, "clust2det.csv")

        run([heliolinc,
             "-imgs",      outim_file,
             "-pairdets",  pairdets_file,
             "-tracklets", trk_file,
             "-trk2det",   trk2det_file,
             "-obspos",    earth_file,
             "-heliodist", hyp_file,
             "-clustrad",  str(CLUSTRAD),
             "-npt",       str(DBSCAN_NPT),
             "-minobsnights", str(MINOBSNIGHTS),
             "-autorun",   "1",
             "-outsum",    outsum_file,
             "-clust2det", c2d_file],
            label="heliolinc", verbose=args.verbose)

        n_raw_clusters = count_data_lines(outsum_file)
        print(f"\nheliolinc: found {n_raw_clusters} raw cluster candidates "
              f"(need ≥ {N_ASTEROIDS})")
        if n_raw_clusters < N_ASTEROIDS:
            print(f"FAIL: heliolinc found only {n_raw_clusters} clusters, "
                  f"expected at least {N_ASTEROIDS}")
            sys.exit(1)

        # Check recovery via clust2det: build a map from paireddet index
        # to asteroid id (using the idstring column), then for each cluster
        # find the dominant asteroid.  A cluster "recovers" an asteroid when
        # >= MINOBSNIGHTS*2 of its detections come from that asteroid.
        from collections import Counter

        # pairdets: header starts with '#'; data rows are 0-indexed.
        det_to_ast = {}
        data_idx = 0
        with open(pairdets_file) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                cols = stripped.split(",")
                if len(cols) > 10:
                    det_to_ast[data_idx] = cols[10]  # idstring column
                data_idx += 1

        # clust2det: cluster_num → Counter of asteroid ids
        cluster_ast: dict = {}
        with open(c2d_file) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split(",")
                if len(parts) < 2:
                    continue
                try:
                    cnum = int(parts[0])
                    dnum = int(parts[1])
                except ValueError:
                    continue
                ast_id = det_to_ast.get(dnum, "unknown")
                cluster_ast.setdefault(cnum, Counter())[ast_id] += 1

        recovered_asteroids = set()
        for counter in cluster_ast.values():
            dominant_ast, count = counter.most_common(1)[0]
            # Require contributions from at least MINOBSNIGHTS nights, i.e.
            # at least MINOBSNIGHTS*2 paired detections from the same object.
            if count >= MINOBSNIGHTS * 2:
                recovered_asteroids.add(dominant_ast)

        print(f"heliolinc recovered objects: {sorted(recovered_asteroids)}")
        if len(recovered_asteroids) < N_ASTEROIDS:
            missing = {a["id"] for a in ASTEROIDS} - recovered_asteroids
            print(f"FAIL: the following asteroids were NOT recovered by "
                  f"heliolinc: {sorted(missing)}")
            sys.exit(1)

        # ------------------------------------------------------------------ #
        # 4. link_purify
        # ------------------------------------------------------------------ #
        lflist_file   = os.path.join(tmpdir, "lflist.txt")
        lp_sum_file   = os.path.join(tmpdir, "LP_clusters.csv")
        lp_c2d_file   = os.path.join(tmpdir, "LP_clust2det.csv")

        # lflist format: one line per heliolinc run → "sumfile clust2detfile"
        with open(lflist_file, "w") as f:
            f.write(f"{outsum_file} {c2d_file}\n")

        run([link_pur,
             "-imgs",    outim_file,
             "-pairdet", pairdets_file,
             "-lflist",  lflist_file,
             "-max_astrom_rms", "2.0",   # generous for synthetic data
             "-outsum",  lp_sum_file,
             "-clust2det", lp_c2d_file],
            label="link_purify", verbose=args.verbose)

        n_final = count_data_lines(lp_sum_file)
        print(f"\nlink_purify: {n_final} final linkages "
              f"(expected {N_ASTEROIDS})")
        if n_final < N_ASTEROIDS:
            print(f"FAIL: link_purify recovered only {n_final} linkages, "
                  f"expected {N_ASTEROIDS}")
            sys.exit(1)

        # ------------------------------------------------------------------ #
        # 5. Result
        # ------------------------------------------------------------------ #
        print(f"\n{'='*60}")
        print(f"PASS — all {N_ASTEROIDS} simulated asteroids recovered.")
        print(f"  make_tracklets : {n_tracklets} tracklets")
        print(f"  heliolinc      : {n_raw_clusters} raw clusters")
        print(f"  link_purify    : {n_final} final linkages")
        print(f"{'='*60}")

    finally:
        if not args.keep_tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            print(f"\nTest files left in: {tmpdir}")


if __name__ == "__main__":
    main()
