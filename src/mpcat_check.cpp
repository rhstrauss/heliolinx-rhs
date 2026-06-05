// mpcat_check.cpp  (June 04, 2026)
//
// Determine exactly which of a set of astrometric measurements already exist
// in the Minor Planet Center observation archive. The heliolinx-native,
// git-installable C++ counterpart of the standalone Python mpcat_check.
//
// It is the rapid, pre-indexed analogue of label_hldet_mpc: instead of loading
// raw MPC 80-column files into RAM (infeasible for the ~40 GB / ~539M-line
// archive), it mmap's the time-sorted binary catalog produced by mpcat_index,
// binary-searches only the MJD window spanned by the input measurements, and
// runs the same 4D (time + sky-position) k-d tree match over that slice.
//
// For each input detection it reports whether a matching MPC observation
// exists (within -matchrad arcsec AND -timerad seconds, optionally requiring
// the same observatory code) and, if so, the matched MPC packed designation.
//
// Usage:
//   mpcat_check -dets in.mpc80 -mpcat mpcat.bin \
//     -matchrad 2.0 -timerad 5.0 -matchobscode 1 -out matched.csv \
//     [-informat mpc80|pairdets] [-colformat colformat.txt] \
//     [-unmatched_out new.csv]
//
// -matchrad is arcseconds (spatial tolerance); -timerad is seconds (temporal
// tolerance). Both must be satisfied for a match. Input format is auto-detected
// from the -dets extension (.mpc80 -> mpc80, else pairdets) unless -informat is
// given.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define DAY_TO_DEG_CONV 24.0
#define COLS_TO_READ 14

static void show_usage()
{
  cerr << "Usage: mpcat_check -dets dets_file -mpcat mpcat.bin -out outfile\n";
  cerr << "       [-matchrad arcsec] [-timerad seconds] [-matchobscode 0/1]\n";
  cerr << "       [-informat mpc80|pairdets] [-colformat colformat.txt]\n";
  cerr << "       [-unmatched_out new.csv]\n";
  cerr << "\n";
  cerr << "Required:\n";
  cerr << "  -dets      Input measurements: MPC 80-column file or heliolinx pairdets CSV\n";
  cerr << "  -mpcat     Pre-indexed binary MPC catalog built by mpcat_index\n";
  cerr << "  -out       Output CSV: one row per input detection with match status\n";
  cerr << "\n";
  cerr << "Optional:\n";
  cerr << "  -matchrad      Spatial match radius, arcsec (default 2.0)\n";
  cerr << "  -timerad       Temporal match tolerance, seconds (default 5.0)\n";
  cerr << "  -matchobscode  Require observatory-code match: 0=no, 1=yes (default 1)\n";
  cerr << "  -informat      Force input format mpc80|pairdets (default: auto by extension)\n";
  cerr << "  -colformat     Column-format file for pairdets input\n";
  cerr << "  -unmatched_out CSV of only the unmatched (genuinely new) input detections\n";
}

int main(int argc, char *argv[])
{
  ifstream instream1;
  ofstream outstream1;
  string colformatfile, stest;
  string dets_file, mpcat_file, outfile, unmatched_outfile, informat;
  vector<hldet> in_dets = {};
  double matchrad = 2.0;   // arcseconds
  double timerad = 5.0;    // seconds
  int matchobscode = 1;
  int verbose = 0;
  int colformatfile_set = 0;
  long i;

  // Column indices (defaults for standard hldet/pairdets format), as in label_hldet_mpc.
  int mjdcol = 1, racol = 2, deccol = 3, magcol = 4;
  int idcol = 11, bandcol = 12, obscodecol = 13;
  int trail_len_col = 5, trail_PA_col = 6, sigmag_col = 7;
  int sig_across_col = 8, sig_along_col = 9;
  int known_obj_col = 14, det_qual_col = 15;

  if(argc < 5) { show_usage(); return(1); }

  i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-dets" || string(argv[i]) == "-d" || string(argv[i]) == "--dets") {
      if(i+1 < argc) { dets_file = argv[++i]; i++; }
      else { cerr << "ERROR: -dets requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-mpcat" || string(argv[i]) == "-cat" || string(argv[i]) == "--mpcat") {
      if(i+1 < argc) { mpcat_file = argv[++i]; i++; }
      else { cerr << "ERROR: -mpcat requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-o" || string(argv[i]) == "--out") {
      if(i+1 < argc) { outfile = argv[++i]; i++; }
      else { cerr << "ERROR: -out requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-matchrad" || string(argv[i]) == "-mrad" || string(argv[i]) == "--matchrad") {
      if(i+1 < argc) { matchrad = stod(argv[++i]); i++; }
      else { cerr << "ERROR: -matchrad requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-timerad" || string(argv[i]) == "-trad" || string(argv[i]) == "--timerad") {
      if(i+1 < argc) { timerad = stod(argv[++i]); i++; }
      else { cerr << "ERROR: -timerad requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-matchobscode" || string(argv[i]) == "-mobscode" || string(argv[i]) == "--matchobscode") {
      if(i+1 < argc) { matchobscode = stoi(argv[++i]); i++; }
      else { cerr << "ERROR: -matchobscode requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-informat" || string(argv[i]) == "-if" || string(argv[i]) == "--informat") {
      if(i+1 < argc) { informat = argv[++i]; i++; }
      else { cerr << "ERROR: -informat requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-cf" || string(argv[i]) == "--colformat") {
      if(i+1 < argc) { colformatfile = argv[++i]; colformatfile_set = 1; i++; }
      else { cerr << "ERROR: -colformat requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-unmatched_out" || string(argv[i]) == "-uout" || string(argv[i]) == "--unmatched_out") {
      if(i+1 < argc) { unmatched_outfile = argv[++i]; i++; }
      else { cerr << "ERROR: -unmatched_out requires an argument\n"; return(1); }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] << "\n";
      i++;
    }
  }

  if(dets_file.empty())  { cerr << "ERROR: no input detections file (-dets)\n"; show_usage(); return(1); }
  if(mpcat_file.empty()) { cerr << "ERROR: no mpcat catalog (-mpcat)\n";       show_usage(); return(1); }
  if(outfile.empty())    { cerr << "ERROR: no output file (-out)\n";           show_usage(); return(1); }

  // Resolve input format: explicit -informat overrides extension auto-detect.
  if(informat.empty()) {
    string lower = dets_file;
    for(auto &ch : lower) ch = tolower(ch);
    if(lower.size() >= 6 && lower.substr(lower.size()-6) == ".mpc80") informat = "mpc80";
    else informat = "pairdets";
  }
  if(informat != "mpc80" && informat != "pairdets") {
    cerr << "ERROR: -informat must be mpc80 or pairdets (got '" << informat << "')\n";
    return(1);
  }

  cout.precision(17);
  cout << "Input detections: " << dets_file << " (format: " << informat << ")\n";
  cout << "MPC catalog: " << mpcat_file << "\n";
  cout << "Spatial match radius: " << matchrad << " arcsec\n";
  cout << "Temporal match tolerance: " << timerad << " seconds\n";
  cout << "Require observatory-code match: " << (matchobscode ? "yes" : "no") << "\n";
  cout << "Output: " << outfile << "\n";

  // ---- Read the input detections ----
  if(informat == "mpc80") {
    int status = read_detection_file_MPC80(dets_file, in_dets);
    if(status != 0) { cerr << "ERROR reading mpc80 file " << dets_file << " (status " << status << ")\n"; return(1); }
  } else {
    // pairdets: optionally read a column-format file (same block as label_hldet_mpc).
    if(colformatfile_set) {
      instream1.open(colformatfile);
      if(!instream1) { cerr << "ERROR: unable to open column format file " << colformatfile << "\n"; return(1); }
      long colreadct = 0;
      while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct < COLS_TO_READ) {
        instream1 >> stest;
        if(stest == "MJDCOL")            { instream1 >> mjdcol;        if(!instream1.fail()) colreadct++; }
        else if(stest == "RACOL")        { instream1 >> racol;         if(!instream1.fail()) colreadct++; }
        else if(stest == "DECCOL")       { instream1 >> deccol;        if(!instream1.fail()) colreadct++; }
        else if(stest == "MAGCOL")       { instream1 >> magcol;        if(!instream1.fail()) colreadct++; }
        else if(stest == "TRAILLENCOL")  { instream1 >> trail_len_col; if(!instream1.fail()) colreadct++; }
        else if(stest == "TRAILPACOL")   { instream1 >> trail_PA_col;  if(!instream1.fail()) colreadct++; }
        else if(stest == "SIGMAGCOL")    { instream1 >> sigmag_col;    if(!instream1.fail()) colreadct++; }
        else if(stest == "SIGACROSSCOL") { instream1 >> sig_across_col;if(!instream1.fail()) colreadct++; }
        else if(stest == "SIGALONGCOL")  { instream1 >> sig_along_col; if(!instream1.fail()) colreadct++; }
        else if(stest == "IDCOL")        { instream1 >> idcol;         if(!instream1.fail()) colreadct++; }
        else if(stest == "BANDCOL")      { instream1 >> bandcol;       if(!instream1.fail()) colreadct++; }
        else if(stest == "OBSCODECOL")   { instream1 >> obscodecol;    if(!instream1.fail()) colreadct++; }
        else if(stest == "KNOWNOBJCOL")  { instream1 >> known_obj_col; if(!instream1.fail()) colreadct++; }
        else if(stest == "DETQUALCOL")   { instream1 >> det_qual_col;  if(!instream1.fail()) colreadct++; }
        else cout << "WARNING: unrecognized string " << stest << " in column format file\n";
      }
      instream1.close();
    }
    int status = read_detection_filemt2(dets_file, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol,
                                        trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col,
                                        known_obj_col, det_qual_col, in_dets, verbose, 1);
    if(status != 0) { cerr << "ERROR reading pairdets file " << dets_file << " (status " << status << ")\n"; return(1); }
  }
  long innum = long(in_dets.size());
  cout << "Read " << innum << " input detections\n";
  if(innum == 0) { cerr << "ERROR: no input detections read\n"; return(1); }

  // ---- mmap the time-sorted binary catalog ----
  // Validate record size against the sidecar header if present.
  {
    ifstream hdrstream(mpcat_file + ".hdr");
    if(hdrstream) {
      string key; long recsize_hdr = -1;
      while(hdrstream >> key) {
        if(key == "recsize") { hdrstream >> recsize_hdr; }
        else { string rest; getline(hdrstream, rest); }
      }
      if(recsize_hdr > 0 && recsize_hdr != long(sizeof(mpcdet))) {
        cerr << "ERROR: mpcat record size mismatch: header says " << recsize_hdr
             << " but this build's sizeof(mpcdet)=" << sizeof(mpcdet)
             << ". The catalog was built with an incompatible mpcdet layout.\n";
        return(1);
      }
    }
  }

  int fd = open(mpcat_file.c_str(), O_RDONLY);
  if(fd < 0) { cerr << "ERROR: can't open mpcat catalog " << mpcat_file << "\n"; return(1); }
  struct stat st;
  if(fstat(fd, &st) != 0) { cerr << "ERROR: fstat failed on " << mpcat_file << "\n"; close(fd); return(1); }
  size_t filesize = size_t(st.st_size);
  if(filesize == 0 || filesize % sizeof(mpcdet) != 0) {
    cerr << "ERROR: " << mpcat_file << " size " << filesize << " is not a whole number of "
         << sizeof(mpcdet) << "-byte records\n";
    close(fd);
    return(1);
  }
  long ncat = long(filesize / sizeof(mpcdet));
  void *map = mmap(nullptr, filesize, PROT_READ, MAP_SHARED, fd, 0);
  if(map == MAP_FAILED) { cerr << "ERROR: mmap failed on " << mpcat_file << "\n"; close(fd); return(1); }
  const mpcdet *cat = reinterpret_cast<const mpcdet*>(map);
  cout << "mmap'd " << ncat << " catalog records (" << double(filesize)/1.0e9 << " GB)\n";
  cout << "Catalog MJD range: " << fixed << setprecision(6) << cat[0].MJD << " to " << cat[ncat-1].MJD << "\n";

  // ---- Binary-search the MJD window spanned by the input detections ----
  double in_tmin = in_dets[0].MJD, in_tmax = in_dets[0].MJD;
  vector<double> mjdvec;
  mjdvec.reserve(innum);
  for(i = 0; i < innum; i++) {
    double m = in_dets[i].MJD;
    mjdvec.push_back(m);
    if(m < in_tmin) in_tmin = m;
    if(m > in_tmax) in_tmax = m;
  }
  double mjdref = dmedian(mjdvec);
  double timerad_days = timerad / SOLARDAY;
  double t0 = in_tmin - timerad_days;
  double t1 = in_tmax + timerad_days;

  // Records are MJD-ascending: lower/upper bound by MJD give the contiguous slice.
  const mpcdet *lo = std::lower_bound(cat, cat + ncat, t0,
                       [](const mpcdet &r, double v){ return r.MJD < v; });
  const mpcdet *hi = std::upper_bound(cat, cat + ncat, t1,
                       [](double v, const mpcdet &r){ return v < r.MJD; });
  long slice_lo = long(lo - cat);
  long slice_hi = long(hi - cat);
  long nslice = slice_hi - slice_lo;
  cout << "Input MJD window [" << fixed << setprecision(6) << t0 << ", " << t1 << "] -> "
       << nslice << " catalog records in slice [" << slice_lo << ", " << slice_hi << ")\n";

  if(nslice == 0) {
    cout << "No catalog records in window; all input detections are NEW.\n";
  }

  // ---- Build a 4D k-d tree over the catalog slice (same geometry as label_hldet_mpc) ----
  double timescale = DAY_TO_DEG_CONV;
  vector<point4d_index> poolvec;
  vector<KD_point4d_index> kdvec;
  if(nslice > 0) {
    poolvec.reserve(nslice);
    for(long j = 0; j < nslice; j++) {
      const mpcdet &r = cat[slice_lo + j];
      poolvec.push_back(point4d_index(
        (r.MJD - mjdref) * timescale / DEGPRAD,
        cos(r.RA / DEGPRAD) * cos(r.Dec / DEGPRAD),
        sin(r.RA / DEGPRAD) * cos(r.Dec / DEGPRAD),
        sin(r.Dec / DEGPRAD),
        j));  // index = slice-local index
    }
    long splitpoint = medind_4d_index(poolvec, 1);
    kdvec.push_back(KD_point4d_index(poolvec[splitpoint], -1, -1, 1, 0));
    kdtree_4d_index(poolvec, 1, splitpoint, 0, kdvec);
    cout << "k-d tree built over slice: " << kdvec.size() << " nodes\n";
  }

  // Combined 4D search radius covering both spatial and temporal tolerances.
  double matchrad_rad = matchrad / ASECPRAD;
  double timerad_rad = timerad_days * timescale / DEGPRAD;
  double search_rad = sqrt(matchrad_rad*matchrad_rad + timerad_rad*timerad_rad);

  // ---- Query each input detection against the slice ----
  vector<int>    in_match(innum, 0);
  vector<long>   in_nmatch(innum, 0);
  vector<string> in_packed(innum, "");
  vector<double> in_sep(innum, -1.0);
  vector<double> in_dt(innum, -1.0);

  long n_in_mpc = 0;
  for(i = 0; i < innum; i++) {
    if(nslice == 0) continue;
    point4d_index q = point4d_index(
      (in_dets[i].MJD - mjdref) * timescale / DEGPRAD,
      cos(in_dets[i].RA / DEGPRAD) * cos(in_dets[i].Dec / DEGPRAD),
      sin(in_dets[i].RA / DEGPRAD) * cos(in_dets[i].Dec / DEGPRAD),
      sin(in_dets[i].Dec / DEGPRAD),
      i);
    vector<long> indexvec;
    kdrange_4d_index(kdvec, q, search_rad, indexvec);

    long nmatch = 0;
    double best_sep = 1.0e30, best_dt = 0.0;
    string best_packed = "";
    for(long k = 0; k < long(indexvec.size()); k++) {
      long sidx = kdvec[indexvec[k]].point.index;   // slice-local index
      const mpcdet &r = cat[slice_lo + sidx];

      double dt_sec = fabs(in_dets[i].MJD - r.MJD) * SOLARDAY;
      if(dt_sec > timerad) continue;

      double cos_sep = cos(in_dets[i].RA/DEGPRAD)*cos(in_dets[i].Dec/DEGPRAD) * cos(r.RA/DEGPRAD)*cos(r.Dec/DEGPRAD)
                     + sin(in_dets[i].RA/DEGPRAD)*cos(in_dets[i].Dec/DEGPRAD) * sin(r.RA/DEGPRAD)*cos(r.Dec/DEGPRAD)
                     + sin(in_dets[i].Dec/DEGPRAD) * sin(r.Dec/DEGPRAD);
      if(cos_sep > 1.0) cos_sep = 1.0;
      if(cos_sep < -1.0) cos_sep = -1.0;
      double sep_arcsec = acos(cos_sep) * ASECPRAD;
      if(sep_arcsec > matchrad) continue;

      if(matchobscode && string(in_dets[i].obscode) != string(r.obscode)) continue;

      nmatch++;
      if(sep_arcsec < best_sep) {
        best_sep = sep_arcsec;
        best_dt = (in_dets[i].MJD - r.MJD) * SOLARDAY;
        best_packed = string(r.packed);
      }
    }
    if(nmatch > 0) {
      in_match[i] = 1;
      in_nmatch[i] = nmatch;
      in_packed[i] = best_packed;
      in_sep[i] = best_sep;
      in_dt[i] = best_dt;
      n_in_mpc++;
    }
  }

  // ---- Report ----
  long n_new = innum - n_in_mpc;
  cout << "\nmpcat_check complete:\n";
  cout << "  Input detections:        " << innum << "\n";
  cout << "  Already in MPC:          " << n_in_mpc << "\n";
  cout << "  Genuinely new:           " << n_new << "\n";
  cout << "  Fraction new:            " << fixed << setprecision(2)
       << (innum ? 100.0*double(n_new)/double(innum) : 0.0) << "%\n";

  // Per-detection output.
  outstream1.open(outfile);
  if(!outstream1) { cerr << "ERROR: can't open output file " << outfile << "\n"; munmap(map, filesize); close(fd); return(1); }
  outstream1 << "#MJD,RA,Dec,mag,obscode,idstring,in_mpc,mpc_packed,n_mpc_match,sep_arcsec,dt_sec\n";
  for(i = 0; i < innum; i++) {
    outstream1 << fixed << setprecision(7) << in_dets[i].MJD << "," << in_dets[i].RA << "," << in_dets[i].Dec << ",";
    outstream1 << fixed << setprecision(3) << in_dets[i].mag << ",";
    outstream1 << in_dets[i].obscode << "," << in_dets[i].idstring << ",";
    outstream1 << in_match[i] << "," << in_packed[i] << "," << in_nmatch[i] << ",";
    if(in_match[i]) outstream1 << fixed << setprecision(4) << in_sep[i] << "," << setprecision(4) << in_dt[i] << "\n";
    else            outstream1 << "-1,-1\n";
  }
  outstream1.close();
  cout << "Wrote " << innum << " rows to " << outfile << "\n";

  // Optional: CSV of only the unmatched (genuinely new) input detections.
  if(!unmatched_outfile.empty()) {
    ofstream ustream(unmatched_outfile);
    if(!ustream) { cerr << "WARNING: can't open unmatched output file " << unmatched_outfile << "\n"; }
    else {
      ustream << "#MJD,RA,Dec,mag,obscode,idstring\n";
      for(i = 0; i < innum; i++) {
        if(in_match[i]) continue;
        ustream << fixed << setprecision(7) << in_dets[i].MJD << "," << in_dets[i].RA << "," << in_dets[i].Dec << ",";
        ustream << fixed << setprecision(3) << in_dets[i].mag << ",";
        ustream << in_dets[i].obscode << "," << in_dets[i].idstring << "\n";
      }
      ustream.close();
      cout << "Wrote " << n_new << " unmatched detections to " << unmatched_outfile << "\n";
    }
  }

  munmap(map, filesize);
  close(fd);
  return(0);
}
