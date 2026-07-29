#!/usr/bin/env python3
"""
reverse_design.py -- Analytically invert the heliolinc (HelioLinC3D) clustering
algorithm: given a Keplerian orbit, a reference epoch MJDref, and the clustering
radius, compute the EXACT heliocentric hypothesis (HelioRad, R_dot, R_dubdot)
that heliolinc needs in its grid to recover that object, plus first-order
"capture half-widths" giving the minimal grid spacing that guarantees recovery.

--------------------------------------------------------------------------------
HELIOLINC HYPOTHESIS CONVENTION (established from source; see reverse_design_design.md)
--------------------------------------------------------------------------------
A grid row is:  r(AU)  rdot(AU/day)  mean_accel
mapped to C++ hlradhyp{ HelioRad; R_dot; R_dubdot } (solarsyst_dyn_geo01.h ~L636).

  HelioRad  = heliocentric distance r at MJDref, in AU.
  R_dot     = dr/dt (radial velocity) at MJDref, in AU/day.
  R_dubdot  = NORMALIZED radial acceleration, dimensionless, with the sign
              convention "POSITIVE ACCELERATION IS INWARD TOWARD THE SUN"
              (calc_heliohypmat.cpp L252-254, L277 `accelnorm = accelk/g0`).

The consumer (solarsyst_dyn_geo01.cpp L39558) converts to a physical radial
acceleration via:
    helioacc_kmday2 = R_dubdot * ( -GMSUN_KM3_SEC2 * SOLARDAY^2 / (r_km)^2 )
                    = R_dubdot * ( -g0 )                 [g0 = GMsun/r^2 > 0]
i.e. the physical outward radial acceleration is  d2r/dt2 = -R_dubdot * g0.
Therefore, for an object whose true outward radial acceleration is `rddot`:
    R_dubdot = -rddot / g0 ,      g0 = mu / r^2   (mu = GMsun).
For a body momentarily on a purely radial free-fall (no transverse motion),
rddot = -g0  =>  R_dubdot = +1  (consistent with the rdot~0 reference rows of
grid_master.txt which sit near mean_accel = 1.0).

heliolinc then uses the quadratic r(t) = r + rdot*dt + 0.5*helioacc*dt^2
(trk2statevec, L33415) to project each tracklet to a 3-D heliocentric position,
forms a velocity, Kepler-integrates the midpoint state to MJDref, and clusters
the resulting 6-D state vectors (position km; velocity km/s * chartimescale)
with a Euclidean radius `clustrad` km (at REF_GEODIST=1 AU; scaled linearly with
geocentric distance).  clustrad is a full 6-D L2 distance.

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
  reverse_design -orbits <orbits.txt> -earth <earth.csv> -mjdref <MJD>
                 [-clustrad 2.0e5] [-baseline_days T] -out <csv>
      -> per-object predictions + analytic capture half-widths to <out>.

  reverse_design -orbits <orbits.txt> -earth <earth.csv> -mjdref <MJD>
                 -emit_grid <obj_index> <gridfile>
      -> writes a single-row hypothesis file for one object (validation probe).
"""
import sys, math, argparse, bisect

# Physical constants (match heliolinx solarsyst_dyn_geo01.h exactly)
K_GAUSS = 0.01720209895               # AU^(3/2)/day; mu = k^2 (AU^3/day^2)
MU      = K_GAUSS * K_GAUSS
AU_KM   = 1.495978700e8               # km   (== AU_KM in the C++)
DAY_S   = 86400.0                     # s    (== SOLARDAY)
TIMECONVSCALE = 4.0                   # heliolinc characteristic-timescale divisor


# ----------------------------------------------------------------------------
# Earth ephemeris (Horizons CSV: JD_TDB, ..., X,Y,Z (km), VX,VY,VZ (km/s))
# ----------------------------------------------------------------------------
def load_earth(path):
    jd = []; cols = [[], [], [], [], [], []]
    insoe = False
    for line in open(path):
        s = line.strip()
        if s.startswith('$$SOE'): insoe = True; continue
        if s.startswith('$$EOE'): break
        if not insoe: continue
        f = [t.strip() for t in s.split(',')]
        try:
            jd.append(float(f[0]))
            for k in range(6):
                cols[k].append(float(f[2 + k]))
        except (ValueError, IndexError):
            continue
    return jd, cols


def earth_state(jd, cols, mjd):
    """Earth heliocentric ecliptic state at MJD -> (pos AU, vel AU/day)."""
    jdq = mjd + 2400000.5
    i = bisect.bisect_left(jd, jdq)
    i = max(1, min(i, len(jd) - 1))
    t0, t1 = jd[i - 1], jd[i]
    w = (jdq - t0) / (t1 - t0)
    out = [cols[k][i - 1] * (1 - w) + cols[k][i] * w for k in range(6)]
    px, py, pz, vx, vy, vz = out
    return ((px / AU_KM, py / AU_KM, pz / AU_KM),
            (vx * DAY_S / AU_KM, vy * DAY_S / AU_KM, vz * DAY_S / AU_KM))


# ----------------------------------------------------------------------------
# Kepler: elements -> heliocentric ecliptic Cartesian state (R, V)
# ----------------------------------------------------------------------------
def kepler_E(M, e, tol=1e-14, itmax=200):
    M = (M + math.pi) % (2 * math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(itmax):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def elements_to_rv(a, e, inc, Om, om, M, mu=MU):
    """Keplerian elements (AU, deg, deg-mean-anomaly) at their own epoch ->
    heliocentric ecliptic (R AU, V AU/day).  Pass M already advanced to the
    epoch you want the state at."""
    d = math.pi / 180.0
    inc *= d; Om *= d; om *= d; M *= d
    E = kepler_E(M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    r = a * (1 - e * cosE)
    # Perifocal position and velocity
    xp = a * (cosE - e)
    yp = a * math.sqrt(1 - e * e) * sinE
    n = math.sqrt(mu / (a ** 3))            # mean motion (rad/day)
    Edot = n / (1 - e * cosE)
    vxp = -a * sinE * Edot
    vyp = a * math.sqrt(1 - e * e) * cosE * Edot
    # Rotation perifocal -> ecliptic
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(om), math.sin(om)
    R11 = cO * cw - sO * sw * ci; R12 = -cO * sw - sO * cw * ci
    R21 = sO * cw + cO * sw * ci; R22 = -sO * sw + cO * cw * ci
    R31 = sw * si;                R32 = cw * si
    Rx = R11 * xp + R12 * yp
    Ry = R21 * xp + R22 * yp
    Rz = R31 * xp + R32 * yp
    Vx = R11 * vxp + R12 * vyp
    Vy = R21 * vxp + R22 * vyp
    Vz = R31 * vxp + R32 * vyp
    return (Rx, Ry, Rz), (Vx, Vy, Vz)


def advance_M(a, M_deg, dt_days):
    """Advance mean anomaly (deg) by dt days for semi-major axis a (AU)."""
    n_deg = math.degrees(math.sqrt(MU / a ** 3))   # deg/day
    return M_deg + n_deg * dt_days


def vdot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def vnorm(a): return math.sqrt(vdot(a, a))


# ----------------------------------------------------------------------------
# Core inversion: orbit -> exact heliolinc hypothesis at MJDref
# ----------------------------------------------------------------------------
def orbit_to_hypothesis(orb, mjdref):
    """orb = dict(a,e,inc,node,argperi,M,epoch).  Returns
    (HelioRad AU, R_dot AU/day, R_dubdot dimensionless, extras)."""
    M_ref = advance_M(orb['a'], orb['M'], mjdref - orb['epoch'])
    R, V = elements_to_rv(orb['a'], orb['e'], orb['inc'],
                          orb['node'], orb['argperi'], M_ref)
    r = vnorm(R)
    v = vnorm(V)
    rdot = vdot(R, V) / r                                   # AU/day (radial velocity)
    # Radial acceleration:  r_ddot = (v^2 - mu/r - rdot^2) / r   (AU/day^2)
    # (from r_ddot = (|V|^2 + R.A - rdot^2)/r with R.A = -mu/r)
    rddot = (v * v - MU / r - rdot * rdot) / r
    g0 = MU / (r * r)                                       # AU/day^2 (>0)
    R_dubdot = -rddot / g0                                  # heliolinc normalized accel
    extras = dict(R=R, V=V, r=r, v=v, rdot=rdot, rddot=rddot, g0=g0)
    return r, rdot, R_dubdot, extras


# ----------------------------------------------------------------------------
# Capture half-widths (first-order).
#
# heliolinc projects each tracklet to a heliocentric distance using the
# quadratic  r_h(t) = r + rdot*(t-tref) + 0.5*helioacc*(t-tref)^2,  then to a
# 3-D position along the (fixed) line of sight.  A perturbation of the
# hypothesis by (dr, drdot, dr_dubdot) shifts the assumed heliocentric distance
# at observation time t by
#     d r_h(t) = dr + drdot*dt + 0.5 * d(helioacc) * dt^2          [dt = t-tref]
# where  helioacc = -R_dubdot * g0  =>  d(helioacc) = -g0 * dr_dubdot.
# To first order the projected 3-D heliocentric position moves by ~ d r_h along
# the heliocentric radial direction (the line-of-sight geometry maps a change in
# r_h to a comparable change in the 3-D point; the proportionality is ~1 for the
# near-radial NEO geometry and we use it as the leading-order estimate).
#
# Two tracklets observed at times t_a, t_b then separate in projected position
# by  |d r_h(t_a) - d r_h(t_b)|.  For the cluster to survive, the spread of the
# 6-D state vectors must stay within clustrad.  Using the half-baseline
# T = 0.5*(MJDmax-MJDmin) as the characteristic |dt|, the per-axis half-width is
# the hypothesis offset that moves the EDGE tracklet by clustrad relative to the
# reference-epoch tracklet:
#     dr           :  d r_h = dr                      -> half-width = clustrad_AU
#     drdot        :  d r_h = drdot * T               -> half-width = clustrad_AU / T
#     dr_dubdot    :  d r_h = 0.5 * g0 * dr_dubdot*T^2 -> half-width = 2*clustrad_AU/(g0*T^2)
# (clustrad_AU = clustrad_km / AU_KM).  These are the positional channel; the
# velocity channel of the 6-D metric gives comparable or looser limits for the
# short NEO arcs, so the position channel is the binding constraint and yields a
# conservative (tighter) half-width -> safe minimal grid spacing = 2*half-width.
# ----------------------------------------------------------------------------
def capture_halfwidths(g0, clustrad_km, baseline_days):
    clustrad_au = clustrad_km / AU_KM
    T = 0.5 * baseline_days                 # half-baseline (days)
    dr = clustrad_au                        # AU
    drdot = clustrad_au / T                 # AU/day
    dr_dubdot = 2.0 * clustrad_au / (g0 * T * T)  # dimensionless
    return dr, drdot, dr_dubdot


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def load_orbits(path):
    """desig a e i node argperi M epoch H G  (header lines start with #)."""
    orbits = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        f = s.split()
        orbits.append(dict(desig=f[0], a=float(f[1]), e=float(f[2]),
                           inc=float(f[3]), node=float(f[4]),
                           argperi=float(f[5]), M=float(f[6]),
                           epoch=float(f[7]),
                           H=float(f[8]) if len(f) > 8 else 99.0,
                           G=float(f[9]) if len(f) > 9 else 0.15))
    return orbits


def outim_baseline(path):
    """Return (MJDmin, MJDmax, MJDref) from an outim file (col1 = MJD)."""
    mjds = []
    for line in open(path):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        try:
            mjds.append(float(s.split()[0]))
        except (ValueError, IndexError):
            continue
    lo, hi = min(mjds), max(mjds)
    return lo, hi, 0.5 * lo + 0.5 * hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-orbits', required=True)
    ap.add_argument('-earth', required=True, help='(loaded for parity; '
                    'inversion needs no Earth state, kept for interface symmetry)')
    ap.add_argument('-mjdref', type=float, required=True)
    ap.add_argument('-clustrad', type=float, default=2.0e5, help='km (6-D L2)')
    ap.add_argument('-baseline_days', type=float, default=None,
                    help='full MJD span of the cadence; if omitted, half-widths '
                         'use a 1-day baseline placeholder')
    ap.add_argument('-out', default=None)
    ap.add_argument('-emit_grid', nargs=2, metavar=('OBJ_INDEX', 'GRIDFILE'),
                    default=None, help='write 1-row hypothesis for one orbit')
    args = ap.parse_args()

    orbits = load_orbits(args.orbits)
    baseline = args.baseline_days if args.baseline_days else 1.0

    # ---- single-object grid-emit mode (validation probe) ----
    if args.emit_grid is not None:
        idx = int(args.emit_grid[0]); gridfile = args.emit_grid[1]
        orb = orbits[idx]
        r, rdot, rdd, ex = orbit_to_hypothesis(orb, args.mjdref)
        with open(gridfile, 'w') as g:
            g.write("#r(AU) rdot(AU/day) mean_accel\n")
            g.write("%.10f %.10f %.10f\n" % (r, rdot, rdd))
        sys.stderr.write("emit_grid obj %d (%s): r=%.6f rdot=%.8f mean_accel=%.6f -> %s\n"
                         % (idx, orb['desig'], r, rdot, rdd, gridfile))
        return

    # ---- full prediction table ----
    out = open(args.out, 'w') if args.out else sys.stdout
    out.write("#obj_index,desig,HelioRad_AU,R_dot_AUday,R_dubdot,"
              "dr_AU,drdot_AUday,dr_dubdot,"
              "min_grid_dr_AU,min_grid_drdot_AUday,min_grid_dr_dubdot\n")
    for i, orb in enumerate(orbits):
        r, rdot, rdd, ex = orbit_to_hypothesis(orb, args.mjdref)
        dr, drdot, drdd = capture_halfwidths(ex['g0'], args.clustrad, baseline)
        out.write("%d,%s,%.8f,%.10f,%.8f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n"
                  % (i, orb['desig'], r, rdot, rdd,
                     dr, drdot, drdd, 2 * dr, 2 * drdot, 2 * drdd))
    if args.out:
        out.close()
        sys.stderr.write("wrote %d predictions to %s\n" % (len(orbits), args.out))


if __name__ == '__main__':
    main()
