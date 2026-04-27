// heliolinc_estimate: predict memory and runtime for heliolinc_lowmem_omp
//
// Accepts the same CLI flags as heliolinc_lowmem_omp, reads input files
// to count detections/tracklets/hypotheses, and prints estimated peak
// memory and rough wall-clock time. Does NOT run the actual algorithm.
//
// Designed for sizing SLURM allocations before submitting jobs.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <omp.h>
#include <sys/stat.h>

static void show_usage()
{
  cerr << "Usage: heliolinc_estimate [same flags as heliolinc_lowmem_omp]\n";
  cerr << "\nReads input files and prints estimated memory/runtime for heliolinc_lowmem_omp.\n";
  cerr << "Does NOT run the actual linking algorithm.\n\n";
  cerr << "Required flags:\n";
  cerr << "  -imgs imfile\n";
  cerr << "  -pairdets paired_detection_file\n";
  cerr << "  -tracklets tracklet_file\n";
  cerr << "  -trk2det tracklet_to_detection_file\n";
  cerr << "  -obspos observer_position_file\n";
  cerr << "  -heliodist heliocentric_hypothesis_file\n";
  cerr << "\nOptional flags (same defaults as heliolinc_lowmem_omp):\n";
  cerr << "  -max_threads N      Number of OpenMP threads (default: system max)\n";
  cerr << "  -useunivar 0|2      Projection mode (default: 0; 2 = TrackletProjCache)\n";
  cerr << "  -mingeodist AU      Min geocentric distance (default: 0.10)\n";
  cerr << "  -maxgeodist AU      Max geocentric distance (default: 100.0)\n";
  cerr << "  -geologstep factor  Log step for geocentric bins (default: 1.5)\n";
  cerr << "  -safety factor      Memory safety multiplier (default: 1.3)\n";
  cerr << "  -nodemem GB         Node memory limit in GB (default: 1450, Klone cpu-g2)\n";
  cerr << "  -nodecpus N         Node CPU limit (default: 192, Klone cpu-g2)\n";
  cerr << "\nAll other heliolinc flags are accepted and ignored (for CLI compatibility).\n";
}

int main(int argc, char *argv[])
{
  HeliolincConfig config;
  string imfile,pairdetfile,trackletfile,trk2detfile,planetfile,accelfile;
  string sumfile = "sumfile_test.csv";
  string clust2detfile = "clust2detfile_test.csv";
  int max_threads = 0;
  double safety_factor = 1.3;
  double node_mem_gb = 1450.0;  // Klone cpu-g2 default
  int node_cpus = 192;          // Klone cpu-g2 default
  long i=0;

  // ---- CLI parsing: identical to heliolinc_lowmem_omp ----
  i=1;
  while(i<argc) {
    if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) { pairdetfile=argv[++i]; i++; }
      else { cerr << "Input paired detection file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    }  else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) { imfile=argv[++i]; i++; }
      else { cerr << "Image file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-tf" || string(argv[i]) == "-trkfile" || string(argv[i]) == "-trackletfile" || string(argv[i]) == "--trackletfile") {
      if(i+1 < argc) { trackletfile=argv[++i]; i++; }
      else { cerr << "Input tracklet file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) { trk2detfile=argv[++i]; i++; }
      else { cerr << "Input trk2det file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-m" || string(argv[i]) == "-mjd" || string(argv[i]) == "-mjdref" || string(argv[i]) == "-MJD" || string(argv[i]) == "--mjd" || string(argv[i]) == "--MJD" || string(argv[i]) == "--modifiedjulianday" || string(argv[i]) == "--ModifiedJulianDay") {
      if(i+1 < argc) { config.MJDref=stod(argv[++i]); i++; }
      else { cerr << "Reference MJD keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-autorun" || string(argv[i]) == "-auto" || string(argv[i]) == "-automjd" || string(argv[i]) == "-autoMJD" || string(argv[i]) == "--auto_mjd" || string(argv[i]) == "--autoMJD" || string(argv[i]) == "--autorun" || string(argv[i]) == "--automjd") {
      if(i+1 < argc) { config.autorun=stoi(argv[++i]); i++; }
      else { cerr << "Autorun keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-obspos" || string(argv[i]) == "-op" || string(argv[i]) == "-obsvec" || string(argv[i]) == "--observer" || string(argv[i]) == "--observer_position" || string(argv[i]) == "--observer_statevec") {
      if(i+1 < argc) { planetfile=argv[++i]; i++; }
      else { cerr << "Observer position file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-heliodist" || string(argv[i]) == "-hd" || string(argv[i]) == "-heliodva" || string(argv[i]) == "-hdva" || string(argv[i]) == "--heliodistvelacc" || string(argv[i]) == "--heliodva") {
      if(i+1 < argc) { accelfile=argv[++i]; i++; }
      else { cerr << "Heliocentric hypothesis file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-clustrad" || string(argv[i]) == "-cr" || string(argv[i]) == "-crad" || string(argv[i]) == "-cluster" || string(argv[i]) == "--clustrad" || string(argv[i]) == "--clusterradius" || string(argv[i]) == "--clusterrad") {
      if(i+1 < argc) { config.clustrad=stod(argv[++i]); i++; }
      else { cerr << "Clustering radius keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-clustchangerad" || string(argv[i]) == "-clustchangerad") {
      if(i+1 < argc) { config.clustchangerad=stod(argv[++i]); i++; }
      else { cerr << "Transition distance keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-npt" || string(argv[i]) == "-npoints" || string(argv[i]) == "-minpts" || string(argv[i]) == "-np" || string(argv[i]) == "--npt" || string(argv[i]) == "--dbscan_npt" || string(argv[i]) == "--DBSCANnpt") {
      if(i+1 < argc) { config.dbscan_npt=stoi(argv[++i]); i++; }
      else { cerr << "DBSCAN npt keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-minobsnights" || string(argv[i]) == "-mon" || string(argv[i]) == "-minobs" || string(argv[i]) == "-minobsnight" || string(argv[i]) == "--minobsnights" || string(argv[i]) == "--minobsnight" || string(argv[i]) == "-minobsn") {
      if(i+1 < argc) { config.minobsnights=stoi(argv[++i]); i++; }
      else { cerr << "Min observing nights keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-mintimespan" || string(argv[i]) == "-mts" || string(argv[i]) == "-minspan" || string(argv[i]) == "-mintspan" || string(argv[i]) == "--mts" || string(argv[i]) == "--mintimespan" || string(argv[i]) == "--mintspan") {
      if(i+1 < argc) { config.mintimespan=stod(argv[++i]); i++; }
      else { cerr << "Minimum time span keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-mingeodist" || string(argv[i]) == "-mingd" || string(argv[i]) == "-mingeo" || string(argv[i]) == "-mingeod" || string(argv[i]) == "--mingeodist" || string(argv[i]) == "--minimum_geocentr_dist") {
      if(i+1 < argc) { config.mingeodist=stod(argv[++i]); i++; }
      else { cerr << "Min geocentric distance keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-maxgeodist" || string(argv[i]) == "-maxgd" || string(argv[i]) == "-maxgeo" || string(argv[i]) == "-maxgeod" || string(argv[i]) == "--maxgeodist" || string(argv[i]) == "--maximum_geocentr_dist") {
      if(i+1 < argc) { config.maxgeodist=stod(argv[++i]); i++; }
      else { cerr << "Max geocentric distance keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-geologstep" || string(argv[i]) == "-gls" || string(argv[i]) == "-geostep" || string(argv[i]) == "-glogstep" || string(argv[i]) == "--geologstep" || string(argv[i]) == "--geodistlogstep" || string(argv[i]) == "--geodiststep") {
      if(i+1 < argc) { config.geologstep=stod(argv[++i]); i++; }
      else { cerr << "Geocentric distance log step keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-mingeoobs" || string(argv[i]) == "-mgo" || string(argv[i]) == "-minobsdist" || string(argv[i]) == "-mindistobs" || string(argv[i]) == "--min_geocentric_obsdist" || string(argv[i]) == "--min_observation_distance" || string(argv[i]) == "--mingeoobs") {
      if(i+1 < argc) { config.mingeoobs=stod(argv[++i]); i++; }
      else { cerr << "Min geocentric distance at observation keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-minimpactpar" || string(argv[i]) == "-mip" || string(argv[i]) == "-minimp" || string(argv[i]) == "-minimppar" || string(argv[i]) == "--minimum_impact_parameter" || string(argv[i]) == "--minimpactpar" || string(argv[i]) == "--min_impact_par") {
      if(i+1 < argc) { config.minimpactpar=stod(argv[++i]); i++; }
      else { cerr << "Min impact parameter keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-mingeofilter" || string(argv[i]) == "-mgf" || string(argv[i]) == "-min_geodist_filter" || string(argv[i]) == "--mingeofilter") {
      if(i+1 < argc) { config.min_geodist_filter=stod(argv[++i]); i++; }
      else { cerr << "Min geocentric distance filter keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-maxbinsize" || string(argv[i]) == "-mbs" || string(argv[i]) == "-max_bin_size" || string(argv[i]) == "--maxbinsize") {
      if(i+1 < argc) { config.max_statevecs_per_bin=stol(argv[++i]); i++; }
      else { cerr << "Max bin size keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-hypinds" || string(argv[i]) == "--hypinds") {
      if(i+1 < argc) { config.hypinds_file=argv[++i]; i++; }
      else { cerr << "ERROR: -hypinds requires a filename argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-useunivar" || string(argv[i]) == "-use_univar" || string(argv[i]) == "-univar" || string(argv[i]) == "-universalvar") {
      if(i+1 < argc) { config.use_univar=stoi(argv[++i]); i++; }
      else { cerr << "use_univar keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-vinf" || string(argv[i]) == "-maxvinf" || string(argv[i]) == "-max_v_inf") {
      if(i+1 < argc) { config.max_v_inf=stod(argv[++i]); i++; }
      else { cerr << "max_v_inf keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) { config.verbose=stoi(argv[++i]); i++; }
      else { cerr << "Verbosity keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-outsum" || string(argv[i]) == "-sum" || string(argv[i]) == "-rms" || string(argv[i]) == "-outrms" || string(argv[i]) == "-or" || string(argv[i]) == "--outsummaryfile" || string(argv[i]) == "--outsum" || string(argv[i]) == "--sum") {
      if(i+1 < argc) { sumfile=argv[++i]; i++; }
      else { cerr << "Output summary file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-clust2det" || string(argv[i]) == "-c2d" || string(argv[i]) == "-clust2detfile" || string(argv[i]) == "--clust2detfile" || string(argv[i]) == "--clust2det" || string(argv[i]) == "--cluster_to_detection") {
      if(i+1 < argc) { clust2detfile=argv[++i]; i++; }
      else { cerr << "Output clust2det file keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-max_threads" || string(argv[i]) == "-maxthreads" || string(argv[i]) == "--max_threads" || string(argv[i]) == "--maxthreads") {
      if(i+1 < argc) { max_threads=stoi(argv[++i]); i++; }
      else { cerr << "Max threads keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-min_RA" || string(argv[i]) == "-minRA") {
      if(i+1 < argc) { config.min_RA=stod(argv[++i]); i++; }
      else { cerr << "min_RA keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-max_RA" || string(argv[i]) == "-maxRA") {
      if(i+1 < argc) { config.max_RA=stod(argv[++i]); i++; }
      else { cerr << "max_RA keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-min_Dec" || string(argv[i]) == "-minDec") {
      if(i+1 < argc) { config.min_Dec=stod(argv[++i]); i++; }
      else { cerr << "min_Dec keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-max_Dec" || string(argv[i]) == "-maxDec") {
      if(i+1 < argc) { config.max_Dec=stod(argv[++i]); i++; }
      else { cerr << "max_Dec keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-safety" || string(argv[i]) == "--safety" || string(argv[i]) == "-safetyfactor") {
      if(i+1 < argc) { safety_factor=stod(argv[++i]); i++; }
      else { cerr << "Safety factor keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-nodemem" || string(argv[i]) == "--nodemem" || string(argv[i]) == "-node_mem") {
      if(i+1 < argc) { node_mem_gb=stod(argv[++i]); i++; }
      else { cerr << "Node memory keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-nodecpus" || string(argv[i]) == "--nodecpus" || string(argv[i]) == "-node_cpus") {
      if(i+1 < argc) { node_cpus=stoi(argv[++i]); i++; }
      else { cerr << "Node CPUs keyword supplied with no corresponding argument\n"; show_usage(); return(1); }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  // ---- Validate required files ----
  if(imfile.size()<=0 || pairdetfile.size()<=0 || trackletfile.size()<=0 ||
     trk2detfile.size()<=0 || planetfile.size()<=0 || accelfile.size()<=0) {
    cerr << "ERROR: all six input files are required\n";
    show_usage();
    return(1);
  }

  // ---- Count lines in input files (lightweight, no data loaded) ----
  // Count non-header lines in each file to determine N_det, N_trk, etc.
  // This avoids loading gigabytes of data just to count elements.
  auto count_lines = [](const string& filepath, bool has_header) -> long {
    ifstream infile(filepath);
    if(!infile.is_open()) {
      cerr << "ERROR: cannot open " << filepath << "\n";
      return -1;
    }
    long count = 0;
    string line;
    if(has_header) getline(infile, line); // skip header
    while(getline(infile, line)) {
      if(line.size() > 0) count++;
    }
    return count;
  };

  // Count non-comment, non-empty lines (for hypothesis file: 4-column, skip lines with norm<=0)
  auto count_hyp_lines = [](const string& filepath) -> long {
    ifstream infile(filepath);
    if(!infile.is_open()) {
      cerr << "ERROR: cannot open " << filepath << "\n";
      return -1;
    }
    long count = 0;
    string line;
    while(getline(infile, line)) {
      if(line.size() == 0 || line[0] == '#') continue;
      // Parse normalization (3rd column) to check if >0
      double r, rdot, norm;
      istringstream iss(line);
      if(iss >> r >> rdot >> norm) {
        if(norm > 0.0) count++;
      }
    }
    return count;
  };

  // pairdets: CSV with header
  long N_det = count_lines(pairdetfile, true);
  if(N_det < 0) return(1);

  // images: space-delimited with header
  long N_img = count_lines(imfile, true);
  if(N_img < 0) return(1);

  // tracklets: CSV with header
  long N_trk = count_lines(trackletfile, true);
  if(N_trk < 0) return(1);

  // trk2det: space-delimited, no header
  long N_t2d = count_lines(trk2detfile, false);
  if(N_t2d < 0) return(1);

  // Earth ephemeris: CSV with header
  long N_earth = count_lines(planetfile, true);
  if(N_earth < 0) return(1);

  // Hypotheses: count valid (norm>0) lines
  long N_hyp = count_hyp_lines(accelfile);
  if(N_hyp < 0) return(1);

  // If hypinds file specified, count only those hypotheses
  if(config.hypinds_file.size() > 0) {
    long hypcount = count_lines(config.hypinds_file, false);
    if(hypcount > 0) {
      cout << "Using hypinds file: " << config.hypinds_file << " (" << hypcount << " hypotheses selected from " << N_hyp << " total)\n";
      N_hyp = hypcount;
    } else {
      cerr << "WARNING: could not read hypinds file " << config.hypinds_file << ", using all " << N_hyp << " hypotheses\n";
    }
  }

  // ---- Determine thread count ----
  int T;
  if(max_threads > 0) {
    T = max_threads;
  } else {
    T = omp_get_max_threads();
  }

  // ---- Get file sizes for I/O time estimate ----
  struct stat st;
  long total_file_bytes = 0;
  const char* input_files[] = {pairdetfile.c_str(), imfile.c_str(), trackletfile.c_str(),
                                trk2detfile.c_str(), planetfile.c_str(), accelfile.c_str()};
  for(int f=0; f<6; f++) {
    if(stat(input_files[f], &st) == 0) total_file_bytes += st.st_size;
  }

  // ---- Compute struct sizes (compile-time accurate) ----
  size_t sz_hldet = sizeof(hldet);
  size_t sz_hlimage = sizeof(hlimage);
  size_t sz_tracklet = sizeof(tracklet);
  size_t sz_longpair = sizeof(longpair);
  size_t sz_hlradhyp = sizeof(hlradhyp);
  size_t sz_EarthState = sizeof(EarthState);
  size_t sz_TrackletProjCache = sizeof(TrackletProjCache);
  size_t sz_point6ix2 = sizeof(point6ix2);
  size_t sz_point3ix2 = sizeof(point3ix2);
  size_t sz_KD_point3ix2 = sizeof(KD_point3ix2);
  size_t sz_shortclust = sizeof(shortclust);
  size_t sz_uint_pair = sizeof(uint_pair);

  // ---- Memory model ----

  // A. Base data (loaded once, shared across threads read-only)
  double M_base = 0.0;
  M_base += (double)N_det * sz_hldet;            // detvec
  M_base += (double)N_img * sz_hlimage;           // image_log
  M_base += (double)N_trk * sz_tracklet;          // tracklets
  M_base += (double)N_t2d * sz_longpair;          // trk2det
  M_base += (double)N_hyp * sz_hlradhyp;          // radhyp
  M_base += (double)N_earth * sz_EarthState;      // earthpos
  M_base += (double)N_trk * sizeof(double);       // sorted_ang_rates
  if(config.use_univar == 2) {
    M_base += (double)N_trk * sz_TrackletProjCache; // trk_cache
  }

  // B. Per-thread transient (worst case: all tracklets qualify)
  long N_sv = 2L * N_trk;  // Two geometric solutions per tracklet (worst case)

  // Number of geocentric distance bins
  int N_geobins = 0;
  {
    double georadcen = config.mingeodist;
    while(georadcen < config.maxgeodist) {
      N_geobins++;
      georadcen *= config.geologstep;
    }
    if(N_geobins < 1) N_geobins = 1;
  }

  // Worst-case largest bin: assume all state vectors land in one bin
  // More typical: roughly uniform, so N_sv / N_geobins
  long N_bin_worst = N_sv;
  long N_bin_typical = N_sv / N_geobins;

  double M_transient_per_thread = 0.0;
  M_transient_per_thread += (double)N_sv * sz_point6ix2;       // allstatevecs
  M_transient_per_thread += (double)N_sv * sizeof(double);     // sv_geodist
  // KD tree + binstatevecs for one bin at a time (typical case)
  double M_kdtree_typical = (double)N_bin_typical * sz_KD_point3ix2;
  double M_binvecs_typical = (double)N_bin_typical * sz_point3ix2;
  // KD tree + binstatevecs worst case
  double M_kdtree_worst = (double)N_bin_worst * sz_KD_point3ix2;
  double M_binvecs_worst = (double)N_bin_worst * sz_point3ix2;

  double M_transient_typical = M_transient_per_thread + M_kdtree_typical + M_binvecs_typical;
  double M_transient_worst = M_transient_per_thread + M_kdtree_worst + M_binvecs_worst;

  // Per-thread output vectors (small, cleared each hypothesis in streaming mode)
  // Include a modest allocation for shortclust + uint_pair per thread
  double M_output_per_thread = 1.0e6; // ~1 MB buffer estimate per thread

  // C. Optimal thread count for the node
  // Solve: node_mem_bytes / safety_factor >= M_base + T_opt * (M_transient_worst + M_output_per_thread)
  double node_mem_bytes = node_mem_gb * 1.0e9;
  double mem_for_threads = node_mem_bytes / safety_factor - M_base;
  int T_opt;
  if(mem_for_threads <= 0) {
    T_opt = 1; // Base data alone exceeds node memory — flag this
  } else {
    T_opt = (int)(mem_for_threads / (M_transient_worst + M_output_per_thread));
    if(T_opt < 1) T_opt = 1;
    if(T_opt > node_cpus) T_opt = node_cpus;
  }

  // D. Total peak memory (for both user-specified T and optimal T_opt)
  double M_peak_typical = M_base + T * (M_transient_typical + M_output_per_thread);
  double M_peak_worst = M_base + T * (M_transient_worst + M_output_per_thread);
  double M_recommended = M_peak_worst * safety_factor;

  double M_peak_opt = M_base + T_opt * (M_transient_worst + M_output_per_thread);
  double M_recommended_opt = M_peak_opt * safety_factor;

  // ---- Runtime model (rough) ----
  // Empirical constants (initial guesses, to be refined with profiling)
  double read_speed_bytes_per_sec = 200.0e6; // 200 MB/s
  double C_proj_sec = 2.0e-6;   // ~2 microseconds per tracklet projection
  double C_kd_sec = 0.5e-6;     // ~0.5 microseconds per KD-tree point operation (build + query)

  double T_io = (double)total_file_bytes / read_speed_bytes_per_sec;
  double T_per_hyp = (double)N_sv * C_proj_sec + (double)N_sv * log2((double)N_sv) * C_kd_sec;
  double T_total = T_io + ((double)N_hyp * T_per_hyp) / T;
  double T_total_opt = T_io + ((double)N_hyp * T_per_hyp) / T_opt;

  // ---- Helper for human-readable sizes ----
  auto fmt_bytes = [](double bytes) -> string {
    char buf[64];
    if(bytes >= 1.0e12) snprintf(buf, sizeof(buf), "%.1f TB", bytes / 1.0e12);
    else if(bytes >= 1.0e9) snprintf(buf, sizeof(buf), "%.2f GB", bytes / 1.0e9);
    else if(bytes >= 1.0e6) snprintf(buf, sizeof(buf), "%.1f MB", bytes / 1.0e6);
    else if(bytes >= 1.0e3) snprintf(buf, sizeof(buf), "%.1f KB", bytes / 1.0e3);
    else snprintf(buf, sizeof(buf), "%.0f B", bytes);
    return string(buf);
  };

  auto fmt_time = [](double seconds) -> string {
    char buf[64];
    if(seconds >= 3600.0) snprintf(buf, sizeof(buf), "~%.1f hr", seconds / 3600.0);
    else if(seconds >= 60.0) snprintf(buf, sizeof(buf), "~%.0f min", seconds / 60.0);
    else snprintf(buf, sizeof(buf), "~%.1f sec", seconds);
    return string(buf);
  };

  auto fmt_count = [](long n) -> string {
    char buf[64];
    if(n >= 1000000) snprintf(buf, sizeof(buf), "%ld (%s%.1fM%s)", n, "\033[1m", (double)n/1.0e6, "\033[0m");
    else if(n >= 1000) snprintf(buf, sizeof(buf), "%ld (%s%.1fK%s)", n, "\033[1m", (double)n/1.0e3, "\033[0m");
    else snprintf(buf, sizeof(buf), "%ld", n);
    return string(buf);
  };

  // ---- Print report ----
  cout << "\n=== heliolinc_estimate ===\n";
  cout << "Estimates for: heliolinc_lowmem_omp (streaming mode)\n\n";

  cout << "Input counts:\n";
  cout << "  Detections:       " << fmt_count(N_det) << "\n";
  cout << "  Images:           " << fmt_count(N_img) << "\n";
  cout << "  Tracklets:        " << fmt_count(N_trk) << "\n";
  cout << "  Trk2det:          " << fmt_count(N_t2d) << "\n";
  cout << "  Hypotheses:       " << fmt_count(N_hyp) << "\n";
  cout << "  Earth ephemeris:  " << fmt_count(N_earth) << "\n";
  cout << "  Threads (OMP):    " << T << "\n";
  cout << "  use_univar:       " << config.use_univar << "\n";
  cout << "  Geocentric bins:  " << N_geobins << " (" << config.mingeodist << " to " << config.maxgeodist << " AU, step " << config.geologstep << ")\n";
  cout << "  Input file total: " << fmt_bytes(total_file_bytes) << "\n\n";

  cout << "Struct sizes (compiled):\n";
  cout << "  hldet:             " << sz_hldet << " B\n";
  cout << "  hlimage:           " << sz_hlimage << " B\n";
  cout << "  tracklet:          " << sz_tracklet << " B\n";
  cout << "  longpair:          " << sz_longpair << " B\n";
  cout << "  TrackletProjCache: " << sz_TrackletProjCache << " B\n";
  cout << "  point6ix2:         " << sz_point6ix2 << " B\n";
  cout << "  point3ix2:         " << sz_point3ix2 << " B\n";
  cout << "  KD_point3ix2:      " << sz_KD_point3ix2 << " B\n";
  cout << "  shortclust:        " << sz_shortclust << " B\n";
  cout << "  uint_pair:         " << sz_uint_pair << " B\n\n";

  cout << "Memory estimate (heliolinc_lowmem_omp streaming):\n";
  cout << "  Base data (shared):          " << fmt_bytes(M_base) << "\n";
  cout << "    detvec:                    " << fmt_bytes((double)N_det * sz_hldet) << "\n";
  cout << "    image_log:                 " << fmt_bytes((double)N_img * sz_hlimage) << "\n";
  cout << "    tracklets:                 " << fmt_bytes((double)N_trk * sz_tracklet) << "\n";
  cout << "    trk2det:                   " << fmt_bytes((double)N_t2d * sz_longpair) << "\n";
  if(config.use_univar == 2) {
    cout << "    trk_cache (univar=2):      " << fmt_bytes((double)N_trk * sz_TrackletProjCache) << "\n";
  }
  cout << "    other (radhyp,earth,rates): " << fmt_bytes((double)N_hyp * sz_hlradhyp + (double)N_earth * sz_EarthState + (double)N_trk * sizeof(double)) << "\n";
  cout << "  Per-thread transient:        " << fmt_bytes(M_transient_typical) << " typical / " << fmt_bytes(M_transient_worst) << " worst\n";
  cout << "    allstatevecs:              " << fmt_bytes((double)N_sv * sz_point6ix2) << "  (N_sv=" << N_sv << " = 2 x N_trk, worst case)\n";
  cout << "    sv_geodist:                " << fmt_bytes((double)N_sv * sizeof(double)) << "\n";
  cout << "    KD tree (1 bin):           " << fmt_bytes(M_kdtree_typical) << " typical / " << fmt_bytes(M_kdtree_worst) << " worst\n";
  cout << "  With " << T << " threads (as specified):\n";
  cout << "    Typical peak:              " << fmt_bytes(M_peak_typical) << "\n";
  cout << "    Worst-case peak:           " << fmt_bytes(M_peak_worst) << "\n";
  cout << "    With " << safety_factor << "x safety:            " << fmt_bytes(M_recommended) << "\n\n";

  cout << "Node constraints: " << node_mem_gb << " GB memory, " << node_cpus << " CPUs\n\n";

  // Optimal thread recommendation
  cout << "\033[1m";
  cout << ">>> RECOMMENDED: -max_threads " << T_opt << "\n";
  cout << "\033[0m";
  cout << "    Worst-case peak:           " << fmt_bytes(M_peak_opt) << "\n";
  cout << "    With " << safety_factor << "x safety:            " << fmt_bytes(M_recommended_opt) << "  (fits in " << node_mem_gb << " GB node)\n";
  if(M_recommended_opt > node_mem_bytes) {
    cout << "\033[1;31m";
    cout << "    WARNING: even 1 thread exceeds node memory! Base data alone = " << fmt_bytes(M_base) << "\n";
    cout << "    Consider reducing input data (highgrading, time-splitting) or using a larger node.\n";
    cout << "\033[0m";
  }
  cout << "\n";

  // SLURM snippet
  cout << "SLURM snippet:\n";
  long slurm_mem_gb = (long)ceil(M_recommended_opt / 1.0e9);
  if(slurm_mem_gb > (long)node_mem_gb) slurm_mem_gb = (long)node_mem_gb;
  cout << "  #SBATCH --cpus-per-task=" << T_opt << "\n";
  cout << "  #SBATCH --mem=" << slurm_mem_gb << "G\n\n";

  cout << "Runtime estimate (ROUGH — refine with profiling):\n";
  cout << "  File I/O:                    " << fmt_time(T_io) << "\n";
  cout << "  With " << T << " threads:            " << fmt_time(T_total) << "  (" << N_hyp << " hyp / " << T << " threads)\n";
  cout << "  With " << T_opt << " threads (optimal):  " << fmt_time(T_total_opt) << "  (" << N_hyp << " hyp / " << T_opt << " threads)\n\n";

  cout << "Notes:\n";
  cout << "  - Memory uses worst-case N_sv = 2 x N_trk (all tracklets qualify every hypothesis)\n";
  cout << "  - Runtime constants are initial guesses; refine with -verbose timing from real runs\n";
  cout << "  - 'Typical' assumes state vectors distribute evenly across geo bins\n";
  cout << "  - 'Worst-case' assumes all state vectors land in one geo bin\n";
  cout << "  - Optimal threads: floor((node_mem/" << safety_factor << " - base) / per_thread_worst), capped at " << node_cpus << " CPUs\n";
  cout << "  - Override node specs with -nodemem GB and -nodecpus N\n";

  return(0);
}
