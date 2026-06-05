// inject_fakes: May 28, 2026
//
// Generate a "pure" fakes-only synthetic detection catalog for the heliolinx
// pipeline: replay a real observing window's cadence and geometry, inject
// synthetic solar-system objects with KNOWN orbits, and emit ONLY the injected
// detections (no real sources, no noise, no stationary sources). The output is
// a canonical 14-column detection catalog that feeds make_tracklets directly,
// plus a truth side-file for recovery scoring.
//
// The tool reuses the heliolinx ephemeris stack (the same routines heliolinc
// uses internally), so injected positions are self-consistent with the linker
// by construction:
//   - Kepler2dyn()         orbital elements -> heliocentric ecliptic state (AU)
//   - read_horizons_csv()  Earth heliocentric ephemeris (km, km/s)
//   - load_image_table()   group a real catalog into images + observer state
//   - read_image_file()    read an existing make_tracklets "outim" manifest
//   - observer state X,Y,Z is heliocentric, in km (as written to the outim file)
//   - statevec_to_celestial()  ecliptic Cartesian difference -> equatorial RA/Dec
//
// UNITS NOTE (critical): Kepler2dyn works in AU / AU/day; the observer state
// (outim X,Y,Z and observer_baryvel01 output) is in km / km/s. We work in AU for
// the asteroid and convert the observer position km->AU before differencing.
// Light-time uses CLIGHT_AUDAY (speed of light in AU/day).
//
// TIME NOTE: observer_baryvel01 (called inside load_image_table) and Kepler2dyn
// expect the observed UTC/UT1 MJD directly; the TT conversion for the ephemeris
// is handled internally. We never pre-add TTDELTAT.

#include "solarsyst_dyn_geo01.h"

#define DEFAULT_IMRAD 3.9        // FOV radius from boresight (deg); ATLAS-ish
#define DEFAULT_FWHM 2.0         // seeing FWHM (arcsec)
#define DEFAULT_M5PCT 95.0       // percentile of per-image mags used as m5
#define DEFAULT_ROLLWIDTH 0.15   // detectability rollover width (mag)
#define DEFAULT_BAND "o"
#define DEFAULT_EXPTIME 30.0     // seconds, if not available from manifest
#define DEFAULT_G 0.15
#define TRAIL_A 0.67             // Veres+2012 detection-loss coefficients
#define TRAIL_B 1.16
#define SNR_AT_M5 5.0            // limiting mag is by definition a 5-sigma cut
#define LIGHTTIME_ITER 3
#define RATE_DT_SEC 60.0         // bracket interval for sky-rate / trail (s)

static void show_usage()
{
  cerr << "USAGE: inject_fakes -earth EarthEphem -obscode ObsCodes\n";
  cerr << "         (-incat real_catalog [-colformat colformatfile] | -outim outim_manifest)\n";
  cerr << "         (-orbits orbitfile | -grid ...) [-population model file]\n";
  cerr << "         -outdets fakes.csv [-outcolformat cf.txt] -truth truth.csv\n\n";
  cerr << "Exposure manifest (need at least one of):\n";
  cerr << "  -incat/-indets  real detection catalog (gives cadence, boresight, per-image m5, band)\n";
  cerr << "  -colformat      column-format file for -incat (default: canonical 14-col)\n";
  cerr << "  -outim/-inimgs  make_tracklets outim file (geometry + observer state + exptime)\n";
  cerr << "                  (m5 unavailable from outim alone: use -alldetections or -maglimit)\n\n";
  cerr << "Orbit source (need at least one of):\n";
  cerr << "  -orbits         orbit list, one per line: desig a e i node argperi M epoch H [G]\n";
  cerr << "                  (a in AU; angles in deg; epoch in MJD)\n";
  cerr << "  -grid           controlled grid: -grid a amin amax na e emin emax ne i imin imax ni H Hmin Hmax nH\n";
  cerr << "  -gridangles     node argperi M (deg) for grid objects (default random per object)\n";
  cerr << "  -gridepoch      epoch MJD for grid orbits (default: manifest midpoint)\n\n";
  cerr << "Detectability / physics:\n";
  cerr << "  -imrad/-fov     FOV radius from boresight, deg (default 3.9)\n";
  cerr << "  -fwhm           seeing FWHM, arcsec (default 2.0)\n";
  cerr << "  -m5pct          percentile of per-image mag used as m5 (default 95)\n";
  cerr << "  -maglimit       fixed limiting mag when no per-image m5 is available\n";
  cerr << "  -rollwidth      detectability rollover width, mag (default 0.15)\n";
  cerr << "  -alldetections  inject every geometrically visible detection (bypass photometric cut)\n";
  cerr << "  -trailed        yes|no|auto (default auto: trail if length > FWHM)\n";
  cerr << "  -exptime        default exposure time (s) if manifest lacks it (default 30)\n";
  cerr << "  -bandcolor      additive V->band color term (default 0)\n";
  cerr << "  -defaultband    band used when none derivable (default o)\n\n";
  cerr << "Output / control:\n";
  cerr << "  -outdets        output fakes-only catalog (canonical 14-col)\n";
  cerr << "  -outcolformat   write the matching colformat file here (optional)\n";
  cerr << "  -truth          truth side-file (orbit <-> detections)\n";
  cerr << "  -labeled        yes|no (default yes: known_obj = object index)\n";
  cerr << "  -nobj           cap the number of injected objects\n";
  cerr << "  -detqual        det_qual value to write (default -1 = unknown)\n";
  cerr << "  -seed           RNG seed (default: time-based)\n";
  cerr << "  -verbose        verbosity\n";
}

// ---------------------------------------------------------------------------
// Local physics / parsing helpers (kept local to avoid touching the shared,
// actively-developed solarsyst_dyn_geo01.cpp on this branch).
// ---------------------------------------------------------------------------

// Mean daily motion (deg/day) for a heliocentric orbit of semimajor axis a (AU).
// KCONST is the Gaussian gravitational constant: mean motion (rad/day) at a=1 AU.
static double mean_daily_motion_degday(double a_au)
{
  return KCONST * pow(a_au, -1.5) * DEGPRAD;
}

// IAU H-G apparent magnitude. r, delta in AU; phase angle in radians.
static double apparent_mag_HG(double H, double G, double r_au, double delta_au, double phase_rad)
{
  double tanhalf = tan(0.5*phase_rad);
  if(tanhalf < 0.0) tanhalf = 0.0;
  double phi1 = exp(-3.33*pow(tanhalf, 0.63));
  double phi2 = exp(-1.87*pow(tanhalf, 1.22));
  double bracket = (1.0-G)*phi1 + G*phi2;
  if(bracket <= 0.0) bracket = 1.0e-30;
  return H + 5.0*log10(r_au*delta_au) - 2.5*log10(bracket);
}

// Trailing detection loss (mag). mu in arcsec/sec, exptime in sec, fwhm in arcsec.
// dm = 1.25*log10(1 + a*x^2/(1+b*x)),  x = mu*t_exp/FWHM   (Veres et al. 2012).
static double trailing_loss_mag(double mu_arcsec_per_sec, double exptime_s, double fwhm_arcsec)
{
  if(fwhm_arcsec <= 0.0 || mu_arcsec_per_sec <= 0.0 || exptime_s <= 0.0) return 0.0;
  double x = mu_arcsec_per_sec*exptime_s/fwhm_arcsec;
  double arg = 1.0 + TRAIL_A*x*x/(1.0 + TRAIL_B*x);
  if(arg <= 0.0) return 0.0;
  return 1.25*log10(arg);
}

// Logistic detection probability in (m_eff - m5). width<=0 -> hard step.
static double detection_probability(double mag_eff, double m5, double width)
{
  if(width <= 0.0) return (mag_eff <= m5) ? 1.0 : 0.0;
  double z = (mag_eff - m5)/width;
  if(z > 40.0) return 0.0;
  if(z < -40.0) return 1.0;
  return 1.0/(1.0 + exp(z));
}

// Column-format descriptor, defaulting to the canonical 14-column layout.
struct ColFmt {
  int idcol, mjdcol, racol, deccol, magcol, bandcol, obscodecol;
  int trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col;
  ColFmt() : idcol(1), mjdcol(2), racol(3), deccol(4), magcol(5), bandcol(6), obscodecol(7),
             trail_len_col(8), trail_PA_col(9), sigmag_col(10), sig_across_col(11),
             sig_along_col(12), known_obj_col(13), det_qual_col(14) {}
};

// Parse a colformat file (same token vocabulary as make_tracklets).
static int parse_colformat_file(const string &fn, ColFmt &cf)
{
  ifstream instream(fn);
  if(!instream) { cerr << "ERROR: cannot open colformat file " << fn << "\n"; return 1; }
  string tok;
  while(instream >> tok) {
    if(tok == "IDCOL") instream >> cf.idcol;
    else if(tok == "MJDCOL") instream >> cf.mjdcol;
    else if(tok == "RACOL") instream >> cf.racol;
    else if(tok == "DECCOL") instream >> cf.deccol;
    else if(tok == "MAGCOL") instream >> cf.magcol;
    else if(tok == "BANDCOL") instream >> cf.bandcol;
    else if(tok == "OBSCODECOL") instream >> cf.obscodecol;
    else if(tok == "TRAILLENCOL") instream >> cf.trail_len_col;
    else if(tok == "TRAILPACOL") instream >> cf.trail_PA_col;
    else if(tok == "SIGMAGCOL") instream >> cf.sigmag_col;
    else if(tok == "SIGACROSSCOL") instream >> cf.sig_across_col;
    else if(tok == "SIGALONGCOL") instream >> cf.sig_along_col;
    else if(tok == "KNOWNOBJCOL") instream >> cf.known_obj_col;
    else if(tok == "DETQUALCOL") instream >> cf.det_qual_col;
    else cerr << "WARNING: unrecognized token '" << tok << "' in colformat file\n";
  }
  return 0;
}

// Read an orbit list: desig a e i node argperi M epoch H [G]. '#' lines skipped.
static int read_orbit_file(const string &fn, vector<asteroid_orbit> &orbits, int verbose)
{
  ifstream instream(fn);
  if(!instream) { cerr << "ERROR: cannot open orbit file " << fn << "\n"; return 1; }
  string line;
  long nread = 0;
  while(getline(instream, line)) {
    if(line.empty() || line[0]=='#') continue;
    // allow comma or whitespace separation
    for(char &c : line) if(c==',') c=' ';
    istringstream ss(line);
    string desig;
    double a,e,i,node,argperi,M,epoch,H,G;
    G = DEFAULT_G;
    if(!(ss >> desig >> a >> e >> i >> node >> argperi >> M >> epoch >> H)) {
      cerr << "WARNING: skipping malformed orbit line: " << line << "\n";
      continue;
    }
    ss >> G; // optional
    if(a <= 0.0 || e < 0.0 || e >= 1.0) {
      cerr << "WARNING: skipping non-elliptical orbit (a=" << a << ", e=" << e << "): " << desig << "\n";
      continue;
    }
    orbits.push_back(asteroid_orbit(desig, a, e, i, node, argperi, M, epoch,
                                    mean_daily_motion_degday(a), H, G));
    nread++;
  }
  if(verbose>0) cout << "Read " << nread << " orbits from " << fn << "\n";
  return 0;
}

// Generate a Cartesian-product grid of orbits over (a, e, i, H).
static void generate_orbit_grid(double amin,double amax,int na,
                                double emin,double emax,int ne,
                                double imin,double imax,int ni,
                                double Hmin,double Hmax,int nH,
                                double epoch, double G,
                                int fixedangles, double node0,double argperi0,double M0,
                                std::mt19937_64 &rng, vector<asteroid_orbit> &orbits)
{
  std::uniform_real_distribution<double> u360(0.0,360.0);
  auto axval = [](double lo,double hi,int n,int k){ return (n<=1)? lo : lo + (hi-lo)*double(k)/double(n-1); };
  long serial=0;
  for(int ia=0; ia<na; ia++) {
    double a = axval(amin,amax,na,ia);
    if(a<=0.0) continue;
    for(int ie=0; ie<ne; ie++) {
      double e = axval(emin,emax,ne,ie);
      if(e<0.0 || e>=1.0) continue;
      for(int ii=0; ii<ni; ii++) {
        double inc = axval(imin,imax,ni,ii);
        for(int ih=0; ih<nH; ih++) {
          double H = axval(Hmin,Hmax,nH,ih);
          double node, argperi, M;
          if(fixedangles) { node=node0; argperi=argperi0; M=M0; }
          else { node=u360(rng); argperi=u360(rng); M=u360(rng); }
          char desig[SHORTSTRINGLEN];
          snprintf(desig, sizeof(desig), "grid%06ld", serial++);
          orbits.push_back(asteroid_orbit(string(desig), a, e, inc, node, argperi, M, epoch,
                                          mean_daily_motion_degday(a), H, G));
        }
      }
    }
  }
}

// Percentile of a set of values (linear, on a sorted copy). pct in [0,100].
static double percentile(vector<double> v, double pct)
{
  if(v.empty()) return 0.0;
  sort(v.begin(), v.end());
  if(pct<=0.0) return v.front();
  if(pct>=100.0) return v.back();
  double idx = (pct/100.0)*(double(v.size())-1.0);
  long lo = long(floor(idx));
  long hi = long(ceil(idx));
  double frac = idx - double(lo);
  return v[lo]*(1.0-frac) + v[hi]*frac;
}

// Topocentric RA/Dec of an orbit at a given observation time, with iterated
// light-time. obs_au is the observer heliocentric ecliptic position in AU.
// Returns r_helio (AU) and delta (AU) via out-params; RA/Dec in degrees.
//
// t_obs_mjd is the OBSERVED (UTC) MJD. Orbital-element epochs are dynamical
// (TT/TDB), and the observer position (from observer_baryvel01) is referenced
// to the TT instant corresponding to this UTC time, so we propagate the
// asteroid at t_obs + TTDELTAT to keep the asteroid-observer geometry in a
// single (dynamical) time system. Validated against JPL Horizons to <0.02".
static int ephem_at(const asteroid_orbit &orb, double t_obs_mjd, const point3d &obs_au,
                    double &RA, double &Dec, double &r_au, double &delta_au)
{
  double t_obs_tt = t_obs_mjd + TTDELTAT/SOLARDAY;
  double t_emit = t_obs_tt;
  point3d astpos = point3d(0,0,0), astvel = point3d(0,0,0);
  vector<double> topo(3,0.0);
  for(int it=0; it<LIGHTTIME_ITER; it++) {
    if(Kepler2dyn(t_emit, orb, astpos, astvel) != 0) return 1;
    topo[0] = astpos.x - obs_au.x;
    topo[1] = astpos.y - obs_au.y;
    topo[2] = astpos.z - obs_au.z;
    delta_au = nvecabs(topo);
    t_emit = t_obs_tt - delta_au/CLIGHT_AUDAY;
  }
  r_au = vecabs3d(astpos);
  if(statevec_to_celestial(topo, RA, Dec) != 0) return 2;
  return 0;
}

// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  string earthfile, obscodefile, incat, colformatfile, outimfile;
  string orbitfile, popmodel, popfile, outdets, outcolformat, truthfile;
  double imrad = DEFAULT_IMRAD;
  double fwhm = DEFAULT_FWHM;
  double m5pct = DEFAULT_M5PCT;
  double maglimit = -1.0;        // <0 => not set
  double rollwidth = DEFAULT_ROLLWIDTH;
  double exptime_default = DEFAULT_EXPTIME;
  double bandcolor = 0.0;
  string defaultband = DEFAULT_BAND;
  int alldetections = 0;
  string trailedmode = "auto";
  int labeled = 1;
  long nobj_cap = -1;
  long detqual = -1;
  long seed = -1;
  int verbose = 0;
  // grid params
  int use_grid = 0;
  double g_amin=0,g_amax=0,g_emin=0,g_emax=0,g_imin=0,g_imax=0,g_Hmin=0,g_Hmax=0;
  int g_na=1,g_ne=1,g_ni=1,g_nH=1;
  int grid_fixedangles = 0;
  double grid_node=0,grid_argperi=0,grid_M=0;
  double gridepoch = -1.0;       // <0 => use manifest midpoint

  if(argc < 2) { show_usage(); return 1; }

  int i = 1;
  while(i < argc) {
    string a = argv[i];
    if(a=="-earth" || a=="-e") { earthfile=argv[++i]; i++; }
    else if(a=="-obscode" || a=="-obscodes") { obscodefile=argv[++i]; i++; }
    else if(a=="-incat" || a=="-indets" || a=="-dets") { incat=argv[++i]; i++; }
    else if(a=="-colformat" || a=="-cf" || a=="-colfmt") { colformatfile=argv[++i]; i++; }
    else if(a=="-outim" || a=="-inimgs" || a=="-imgs") { outimfile=argv[++i]; i++; }
    else if(a=="-orbits" || a=="-orbit") { orbitfile=argv[++i]; i++; }
    else if(a=="-population" || a=="-pop") { popmodel=argv[++i]; popfile=argv[++i]; i++; }
    else if(a=="-nobj") { nobj_cap=stol(argv[++i]); i++; }
    else if(a=="-imrad" || a=="-fov") { imrad=stod(argv[++i]); i++; }
    else if(a=="-fwhm") { fwhm=stod(argv[++i]); i++; }
    else if(a=="-m5pct") { m5pct=stod(argv[++i]); i++; }
    else if(a=="-maglimit") { maglimit=stod(argv[++i]); i++; }
    else if(a=="-rollwidth") { rollwidth=stod(argv[++i]); i++; }
    else if(a=="-alldetections" || a=="-allvisible") { alldetections=1; i++; }
    else if(a=="-trailed") { trailedmode=argv[++i]; i++; }
    else if(a=="-exptime") { exptime_default=stod(argv[++i]); i++; }
    else if(a=="-bandcolor") { bandcolor=stod(argv[++i]); i++; }
    else if(a=="-defaultband") { defaultband=argv[++i]; i++; }
    else if(a=="-outdets" || a=="-out") { outdets=argv[++i]; i++; }
    else if(a=="-outcolformat") { outcolformat=argv[++i]; i++; }
    else if(a=="-truth") { truthfile=argv[++i]; i++; }
    else if(a=="-labeled") { labeled = (string(argv[++i])=="no"||string(argv[i])=="0")?0:1; i++; }
    else if(a=="-detqual") { detqual=stol(argv[++i]); i++; }
    else if(a=="-seed") { seed=stol(argv[++i]); i++; }
    else if(a=="-verbose" || a=="-v") { verbose=stoi(argv[++i]); i++; }
    else if(a=="-gridangles") { grid_fixedangles=1; grid_node=stod(argv[++i]); grid_argperi=stod(argv[++i]); grid_M=stod(argv[++i]); i++; }
    else if(a=="-gridepoch") { gridepoch=stod(argv[++i]); i++; }
    else if(a=="-grid") {
      use_grid = 1;
      // -grid a amin amax na e emin emax ne i imin imax ni H Hmin Hmax nH
      string ax;
      while(i+1<argc && argv[i+1][0]!='-') {
        ax = argv[++i];
        if(ax=="a") { g_amin=stod(argv[++i]); g_amax=stod(argv[++i]); g_na=stoi(argv[++i]); }
        else if(ax=="e") { g_emin=stod(argv[++i]); g_emax=stod(argv[++i]); g_ne=stoi(argv[++i]); }
        else if(ax=="i") { g_imin=stod(argv[++i]); g_imax=stod(argv[++i]); g_ni=stoi(argv[++i]); }
        else if(ax=="H") { g_Hmin=stod(argv[++i]); g_Hmax=stod(argv[++i]); g_nH=stoi(argv[++i]); }
        else { cerr << "ERROR: unknown -grid axis '" << ax << "'\n"; return 1; }
      }
      i++;
    }
    else { cerr << "ERROR: unrecognized argument '" << a << "'\n"; show_usage(); return 1; }
  }

  // ---- validation ----
  if(earthfile.empty() || obscodefile.empty()) {
    cerr << "ERROR: -earth and -obscode are required\n"; show_usage(); return 1;
  }
  if(incat.empty() && outimfile.empty()) {
    cerr << "ERROR: supply at least one of -incat or -outim (exposure manifest)\n"; return 1;
  }
  if(orbitfile.empty() && !use_grid && popmodel.empty()) {
    cerr << "ERROR: supply at least one orbit source (-orbits, -grid, or -population)\n"; return 1;
  }
  if(outdets.empty()) { cerr << "ERROR: -outdets is required\n"; return 1; }
  if(incat.empty() && !alldetections && maglimit<0.0) {
    cerr << "ERROR: no per-image m5 available from an outim-only manifest; supply -incat,\n"
         << "       or run with -alldetections, or set a fixed -maglimit.\n";
    return 1;
  }

  if(seed < 0) seed = (long)std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937_64 rng((unsigned long long)seed);
  std::uniform_real_distribution<double> uni01(0.0,1.0);
  std::normal_distribution<double> gauss(0.0,1.0);
  cout << "inject_fakes: RNG seed = " << seed << "\n";

  // ---- read observatory codes + Earth ephemeris ----
  vector<observatory> observatory_list;
  if(read_obscode_file2(obscodefile, observatory_list, verbose) != 0) {
    cerr << "ERROR reading obscode file " << obscodefile << "\n"; return 1;
  }
  vector<double> EarthMJD;
  vector<point3d> Earthpos, Earthvel;
  if(read_horizons_csv(earthfile, EarthMJD, Earthpos, Earthvel) != 0) {
    cerr << "ERROR reading Earth ephemeris " << earthfile << "\n"; return 1;
  }
  cout << "Read " << EarthMJD.size() << " Earth ephemeris points and "
       << observatory_list.size() << " observatory codes.\n";

  // ---- build exposure manifest (img_log) + per-image m5 + band ----
  vector<hlimage> img_log;
  vector<hldet> realdets;        // only used on the catalog path (for m5/band)
  vector<double> m5_img;         // per-image limiting mag (catalog path)
  vector<string> band_img;       // per-image band

  if(!outimfile.empty()) {
    if(read_image_file(outimfile, img_log) != 0) {
      cerr << "ERROR reading outim file " << outimfile << "\n"; return 1;
    }
    cout << "Read " << img_log.size() << " images from outim file " << outimfile << "\n";
  }

  if(!incat.empty()) {
    ColFmt cf;
    if(!colformatfile.empty() && parse_colformat_file(colformatfile, cf)!=0) return 1;
    if(read_detection_filemt2(incat, cf.mjdcol, cf.racol, cf.deccol, cf.magcol, cf.idcol,
                              cf.bandcol, cf.obscodecol, cf.trail_len_col, cf.trail_PA_col,
                              cf.sigmag_col, cf.sig_across_col, cf.sig_along_col,
                              cf.known_obj_col, cf.det_qual_col, realdets, verbose, 1) != 0) {
      cerr << "ERROR reading detection catalog " << incat << "\n"; return 1;
    }
    sort(realdets.begin(), realdets.end(), early_hldet());
    cout << "Read " << realdets.size() << " real detections from " << incat << "\n";

    if(outimfile.empty()) {
      // Catalog-only: build the manifest from the catalog exactly as the
      // pipeline does (image grouping + boresight + heliocentric observer state).
      if(load_image_table(img_log, realdets, 0.0, observatory_list, EarthMJD, Earthpos, Earthvel) != 0) {
        cerr << "ERROR: load_image_table failed\n"; return 1;
      }
    } else {
      // outim supplied geometry; match the catalog to those images (by MJD/obscode)
      // only to populate startind/endind for the per-image m5 + band estimate.
      // load_image_indices compares MJD differences directly, so tol is in DAYS.
      load_image_indices(img_log, realdets, IMAGETIMETOL/SOLARDAY, 1);
    }
  }

  if(img_log.empty()) { cerr << "ERROR: no images in manifest\n"; return 1; }

  // Compute the heliocentric observer state (km, km/s) for every image via the
  // same routine the pipeline uses. This is essential for the outim path:
  // read_image_file does NOT parse the X,Y,Z,VX,VY,VZ columns (it leaves them 0),
  // which would otherwise collapse the topocentric vector to a heliocentric one.
  // For the catalog path load_image_table already set these; recomputing is
  // identical and keeps both paths uniform and correct.
  {
    long nbad = 0;
    for(size_t k=0; k<img_log.size(); k++) {
      double obslon, plxcos, plxsin;
      if(obscode_lookup(observatory_list, img_log[k].obscode, obslon, plxcos, plxsin) != 0) {
        nbad++; continue;
      }
      point3d opos(0,0,0), ovel(0,0,0);
      observer_baryvel01(img_log[k].MJD, 5, obslon, plxcos, plxsin,
                         EarthMJD, Earthpos, Earthvel, opos, ovel);
      img_log[k].X = opos.x; img_log[k].Y = opos.y; img_log[k].Z = opos.z;
      img_log[k].VX = ovel.x; img_log[k].VY = ovel.y; img_log[k].VZ = ovel.z;
    }
    if(nbad>0) cerr << "WARNING: obscode_lookup failed for " << nbad << " images (observer state left at 0)\n";
  }

  // exposure time: prefer manifest value; fall back to default
  for(size_t k=0; k<img_log.size(); k++)
    if(!(img_log[k].exptime > 0.0)) img_log[k].exptime = exptime_default;

  // per-image m5 + band from the real catalog (if available)
  m5_img.assign(img_log.size(), maglimit);
  band_img.assign(img_log.size(), defaultband);
  if(!realdets.empty()) {
    for(size_t k=0; k<img_log.size(); k++) {
      long s = img_log[k].startind, e = img_log[k].endind;
      if(s<0 || e<=s || e>(long)realdets.size()) continue;
      vector<double> mags;
      std::unordered_map<string,long> bandcount;
      for(long d=s; d<e; d++) {
        if(isnormal(realdets[d].mag)) mags.push_back(realdets[d].mag);
        bandcount[string(realdets[d].band)]++;
      }
      if(!mags.empty()) m5_img[k] = percentile(mags, m5pct);
      long best=-1; for(auto &p : bandcount) if(p.second>best){best=p.second; band_img[k]=p.first;}
    }
    cout << "Derived per-image m5 (" << m5pct << "th pct) and band from the real catalog.\n";
  }

  // manifest time span (for default grid epoch)
  double mjd_min=img_log.front().MJD, mjd_max=img_log.front().MJD;
  for(auto &im : img_log){ if(im.MJD<mjd_min)mjd_min=im.MJD; if(im.MJD>mjd_max)mjd_max=im.MJD; }
  double mjd_mid = 0.5*(mjd_min+mjd_max);
  cout << "Manifest spans MJD " << fixed << setprecision(5) << mjd_min << " to " << mjd_max
       << " (" << img_log.size() << " images)\n";

  // ---- assemble orbit set ----
  vector<asteroid_orbit> orbits;
  if(!orbitfile.empty() && read_orbit_file(orbitfile, orbits, verbose)!=0) return 1;
  if(use_grid) {
    double ge = (gridepoch>0.0) ? gridepoch : mjd_mid;
    generate_orbit_grid(g_amin,g_amax,g_na, g_emin,g_emax,g_ne, g_imin,g_imax,g_ni,
                        g_Hmin,g_Hmax,g_nH, ge, DEFAULT_G,
                        grid_fixedangles, grid_node,grid_argperi,grid_M, rng, orbits);
    cout << "Generated grid: " << (long)g_na*g_ne*g_ni*g_nH << " orbits (epoch MJD " << ge << ")\n";
  }
  if(!popmodel.empty()) {
    cerr << "WARNING: -population (" << popmodel << ") sampling is not implemented in v1; ignoring.\n";
  }
  if(orbits.empty()) { cerr << "ERROR: no orbits to inject\n"; return 1; }
  if(nobj_cap>=0 && (long)orbits.size()>nobj_cap) orbits.erase(orbits.begin()+nobj_cap, orbits.end());
  cout << "Injecting " << orbits.size() << " synthetic objects.\n";

  int trail_yes = (trailedmode=="yes"||trailedmode=="1");
  int trail_no  = (trailedmode=="no"||trailedmode=="0");

  // ---- main injection loop ----
  vector<hldet> fakedets;
  // truth bookkeeping
  vector<long> tr_nvis(orbits.size(),0), tr_ndet(orbits.size(),0);
  vector< vector<long> > tr_imgs(orbits.size());
  double dt_days = RATE_DT_SEC/SOLARDAY;

  for(size_t objct=0; objct<orbits.size(); objct++) {
    const asteroid_orbit &orb = orbits[objct];
    long serial = 0;
    for(size_t k=0; k<img_log.size(); k++) {
      const hlimage &im = img_log[k];
      point3d obs_au = point3d(im.X/AU_KM, im.Y/AU_KM, im.Z/AU_KM);

      double RA,Dec,r_au,delta_au;
      if(ephem_at(orb, im.MJD, obs_au, RA, Dec, r_au, delta_au)!=0) continue; // unbound/bad

      // FOV / visibility cut
      if(distradec01(RA, Dec, im.RA, im.Dec) > imrad) continue;
      tr_nvis[objct]++;

      // phase angle (Sun-object-observer), s = Sun-observer distance
      double s_au = vecabs3d(obs_au);
      double cosalpha = (r_au*r_au + delta_au*delta_au - s_au*s_au)/(2.0*r_au*delta_au);
      if(cosalpha>1.0) cosalpha=1.0; if(cosalpha<-1.0) cosalpha=-1.0;
      double phase = acos(cosalpha);

      double Vmag = apparent_mag_HG(orb.H, orb.G, r_au, delta_au, phase);
      double obsmag = Vmag + bandcolor;

      // sky rate + trail geometry from a second epoch dt later. The observer is
      // advanced with its velocity (km/s -> AU) so the apparent rate includes
      // Earth's reflex/parallactic motion, which DOMINATES the apparent rate of
      // slow movers (MBAs) near opposition.
      double RA2,Dec2,r2,d2;
      double mu_arcsec_s = 0.0, trail_len = 0.0, trail_PA = 90.0;
      point3d obs_au2 = point3d(obs_au.x + im.VX*RATE_DT_SEC/AU_KM,
                                obs_au.y + im.VY*RATE_DT_SEC/AU_KM,
                                obs_au.z + im.VZ*RATE_DT_SEC/AU_KM);
      if(ephem_at(orb, im.MJD+dt_days, obs_au2, RA2, Dec2, r2, d2)==0) {
        double dist_deg, pa_deg;
        distradec02(RA, Dec, RA2, Dec2, &dist_deg, &pa_deg);
        mu_arcsec_s = dist_deg*3600.0/RATE_DT_SEC;
        trail_len = mu_arcsec_s*im.exptime;            // arcsec
        trail_PA = pa_deg;
        while(trail_PA>=180.0) trail_PA-=180.0;         // trail axis: mod 180
        while(trail_PA<0.0) trail_PA+=180.0;
      }

      // detectability
      double m5 = m5_img[k];
      double dm_trail = trailing_loss_mag(mu_arcsec_s, im.exptime, fwhm);
      double mag_eff = obsmag + dm_trail;
      double SNR;
      if(!alldetections) {
        if(!(m5>0.0)) { /* no m5: treat as always-detected */ }
        else {
          double p = detection_probability(mag_eff, m5, rollwidth);
          if(uni01(rng) > p) continue; // not detected
        }
      }
      if(m5>0.0) SNR = SNR_AT_M5*pow(10.0, 0.4*(m5 - mag_eff));
      else SNR = 100.0; // bright/undefined-m5: assume high SNR
      if(SNR < 1.0) SNR = 1.0;
      if(SNR > 1.0e4) SNR = 1.0e4;

      // trailed decision
      int trailed = trail_yes || (!trail_no && trail_len > fwhm);

      // anisotropic position scatter (arcsec) in (along, cross) then rotate by PA
      double sig_along, sig_cross;
      if(trailed) {
        sig_along = trail_len/(2.0*sqrt(3.0)*SNR);
        sig_cross = fwhm/(2.0*sqrt(2.0)*SNR);
      } else {
        sig_along = sig_cross = fwhm/(2.0*sqrt(2.0)*SNR);
      }
      double da = gauss(rng)*sig_along;
      double dc = gauss(rng)*sig_cross;
      double pa_rad = trail_PA/DEGPRAD;
      double offE = da*sin(pa_rad) + dc*cos(pa_rad);   // arcsec, +East
      double offN = da*cos(pa_rad) - dc*sin(pa_rad);   // arcsec, +North
      double cosdec = cos(Dec/DEGPRAD); if(fabs(cosdec)<1e-8) cosdec=1e-8;
      double RA_out = RA + (offE/3600.0)/cosdec;
      double Dec_out = Dec + offN/3600.0;

      double sigmag = 1.0857/SNR;

      double tl_out = trailed ? trail_len : 0.0;
      double tp_out = trailed ? trail_PA : 90.0;

      char ids[SHORTSTRINGLEN];
      snprintf(ids, sizeof(ids), "fk%ld_%ld", (long)objct, serial++);
      string band = band_img[k];
      if(band.size() >= MINSTRINGLEN) band = band.substr(0, MINSTRINGLEN-1);
      string obsc = string(im.obscode);

      long kobj = labeled ? (long)objct : -1;
      fakedets.push_back(hldet(im.MJD, RA_out, Dec_out, (float)obsmag, (float)tl_out, (float)tp_out,
                               (float)sigmag, (float)sig_cross, (float)sig_along, -1,
                               string(ids), band, obsc, kobj, detqual, (long)fakedets.size()));
      tr_ndet[objct]++;
      tr_imgs[objct].push_back((long)k);
    }
  }

  // sort by MJD and re-index
  sort(fakedets.begin(), fakedets.end(), early_hldet());
  for(size_t k=0; k<fakedets.size(); k++) fakedets[k].index = (long)k;

  // ---- write the fakes-only catalog (canonical 14-col) ----
  {
    FILE *fp = fopen(outdets.c_str(), "w");
    if(!fp) { cerr << "ERROR: cannot open output " << outdets << "\n"; return 1; }
    fprintf(fp, "#ID,MJD,RA,Dec,Mag,Band,ObsCode,trail_len,trail_PA,sigmag,sig_across,sig_along,known_obj,det_qual\n");
    for(const hldet &d : fakedets) {
      fprintf(fp, "%s,%.7f,%.7f,%.7f,%.4f,%s,%s,%.2f,%.2f,%.4f,%.3f,%.3f,%ld,%ld\n",
              d.idstring, d.MJD, d.RA, d.Dec, d.mag, d.band, d.obscode,
              d.trail_len, d.trail_PA, d.sigmag, d.sig_across, d.sig_along,
              d.known_obj, d.det_qual);
    }
    fclose(fp);
    cout << "Wrote " << fakedets.size() << " synthetic detections to " << outdets << "\n";
  }

  if(!outcolformat.empty()) {
    FILE *fp = fopen(outcolformat.c_str(), "w");
    if(fp) {
      fprintf(fp, "IDCOL 1\nMJDCOL 2\nRACOL 3\nDECCOL 4\nMAGCOL 5\nBANDCOL 6\nOBSCODECOL 7\n"
                  "TRAILLENCOL 8\nTRAILPACOL 9\nSIGMAGCOL 10\nSIGACROSSCOL 11\nSIGALONGCOL 12\n"
                  "KNOWNOBJCOL 13\nDETQUALCOL 14\n");
      fclose(fp);
      cout << "Wrote output colformat to " << outcolformat << "\n";
    }
  }

  // ---- write the truth side-file ----
  if(!truthfile.empty()) {
    FILE *fp = fopen(truthfile.c_str(), "w");
    if(!fp) { cerr << "ERROR: cannot open truth file " << truthfile << "\n"; return 1; }
    fprintf(fp, "#obj_index,desig,known_obj,a_au,e,i_deg,node_deg,argperi_deg,M_deg,epoch_mjd,H,G,n_visible,n_detected,image_indices\n");
    for(size_t o=0; o<orbits.size(); o++) {
      const asteroid_orbit &orb = orbits[o];
      fprintf(fp, "%ld,%s,%ld,%.10f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%ld,%ld,",
              (long)o, orb.desig.c_str(), labeled?(long)o:-1,
              orb.semimaj_axis, orb.eccentricity, orb.inclination, orb.long_ascend_node,
              orb.arg_perihelion, orb.mean_anom, orb.mjd_epoch, orb.H, orb.G,
              tr_nvis[o], tr_ndet[o]);
      for(size_t j=0; j<tr_imgs[o].size(); j++) fprintf(fp, "%s%ld", j?";":"", tr_imgs[o][j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    cout << "Wrote truth side-file " << truthfile << "\n";
  }

  // ---- summary ----
  long objs_with_dets=0, total_vis=0, total_det=0;
  for(size_t o=0;o<orbits.size();o++){ total_vis+=tr_nvis[o]; total_det+=tr_ndet[o]; if(tr_ndet[o]>0)objs_with_dets++; }
  cout << "\n=== inject_fakes summary ===\n";
  cout << "objects injected:           " << orbits.size() << "\n";
  cout << "objects with >=1 detection: " << objs_with_dets << "\n";
  cout << "total visible (in FOV):     " << total_vis << "\n";
  cout << "total detected (written):   " << total_det << "\n";
  if(orbits.size()>0)
    cout << "mean detections/object:     " << double(total_det)/double(orbits.size()) << "\n";
  return 0;
}
