// March 2026: merge_det_catalogs.cpp:
// Reads a list file whose lines each contain two fields:
//   catalog_file  colformat_file
// Each catalog_file is a detection catalog in the heliolinc detection format
// (arbitrary column order specified by the accompanying colformat_file).
// All detections from every catalog are merged, sorted by MJD, and written
// to a single output catalog in the canonical 14-column format below, which
// can be fed directly into make_tracklets with a matching colformat file.
//
// Output column order (1-indexed):
//   1  ID
//   2  MJD
//   3  RA         (degrees)
//   4  Dec        (degrees)
//   5  Mag
//   6  Band
//   7  ObsCode
//   8  trail_len  (arcsec; 0 if not in source catalog)
//   9  trail_PA   (degrees; 90 if not in source catalog)
//  10  sigmag     (mag; 9.999 if not in source catalog)
//  11  sig_across (arcsec; 1.0 if not in source catalog)
//  12  sig_along  (arcsec; 1.0 if not in source catalog)
//  13  known_obj  (-1 if not in source catalog)
//  14  det_qual   (-1 if not in source catalog)
//
// Performance notes:
//  - Catalogs are read in parallel (one OpenMP thread per catalog).
//  - Each sub-catalog is sorted independently; if already time-ordered
//    (the common case for telescope data) the sort is skipped after an
//    O(N) is_sorted check.
//  - The sorted sub-catalogs are merged with a min-heap k-way merge:
//    O(N log k) rather than O(N log N) for a global re-sort.
//  - Output is written with fprintf for minimal formatting overhead.
//
// Usage:
//   merge_det_catalogs -catlist catlist_file -out output_file
//   [-verbose verbosity_level]

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <queue>
#include <tuple>

#define IDCOL       1
#define MJDCOL      2
#define RACOL       3
#define DECCOL      4
#define MAGCOL      5
#define BANDCOL     6
#define OBSCODECOL  7
#define COLS_TO_READ 14

static void show_usage()
{
  cerr << "Usage: merge_det_catalogs -catlist catlist_file -out output_file\n";
  cerr << "\nThe catlist_file has one line per input catalog:\n";
  cerr << "  catalog_file  colformat_file\n";
  cerr << "\nThe output is a single time-sorted catalog in the canonical\n";
  cerr << "14-column format (ID,MJD,RA,Dec,Mag,Band,ObsCode,trail_len,\n";
  cerr << "trail_PA,sigmag,sig_across,sig_along,known_obj,det_qual).\n";
  cerr << "Use this colformat with the output and make_tracklets:\n";
  cerr << "  IDCOL 1\n  MJDCOL 2\n  RACOL 3\n  DECCOL 4\n";
  cerr << "  MAGCOL 5\n  BANDCOL 6\n  OBSCODECOL 7\n";
  cerr << "  TRAILLENCOL 8\n  TRAILPACOL 9\n  SIGMAGCOL 10\n";
  cerr << "  SIGACROSSCOL 11\n  SIGALONGCOL 12\n";
  cerr << "  KNOWNOBJCOL 13\n  DETQUALCOL 14\n";
}

// Parse a colformat file into per-column index variables.
// Returns number of recognised keywords found.
static int parse_colformat(const string &colformatfile,
                           int &idcol, int &mjdcol, int &racol, int &deccol,
                           int &magcol, int &bandcol, int &obscodecol,
                           int &trail_len_col, int &trail_PA_col,
                           int &sigmag_col, int &sig_across_col,
                           int &sig_along_col,
                           int &known_obj_col, int &det_qual_col)
{
  ifstream cfstream(colformatfile);
  if(!cfstream) return -1;
  int colreadct = 0;
  string stest;
  while(!cfstream.eof() && !cfstream.fail() && !cfstream.bad()
        && colreadct < COLS_TO_READ) {
    cfstream >> stest;
    if     (stest == "MJDCOL")      { cfstream >> mjdcol;        colreadct++; }
    else if(stest == "RACOL")       { cfstream >> racol;         colreadct++; }
    else if(stest == "DECCOL")      { cfstream >> deccol;        colreadct++; }
    else if(stest == "MAGCOL")      { cfstream >> magcol;        colreadct++; }
    else if(stest == "TRAILLENCOL") { cfstream >> trail_len_col; colreadct++; }
    else if(stest == "TRAILPACOL")  { cfstream >> trail_PA_col;  colreadct++; }
    else if(stest == "SIGMAGCOL")   { cfstream >> sigmag_col;    colreadct++; }
    else if(stest == "SIGACROSSCOL"){ cfstream >> sig_across_col;colreadct++; }
    else if(stest == "SIGALONGCOL") { cfstream >> sig_along_col; colreadct++; }
    else if(stest == "IDCOL")       { cfstream >> idcol;         colreadct++; }
    else if(stest == "BANDCOL")     { cfstream >> bandcol;       colreadct++; }
    else if(stest == "OBSCODECOL")  { cfstream >> obscodecol;    colreadct++; }
    else if(stest == "KNOWNOBJCOL") { cfstream >> known_obj_col; colreadct++; }
    else if(stest == "DETQUALCOL")  { cfstream >> det_qual_col;  colreadct++; }
  }
  return colreadct;
}

int main(int argc, char *argv[])
{
  string catlistfile, outfile;
  int verbose = 0;

  if(argc < 5) {
    show_usage();
    return(1);
  }

  int i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-catlist" || string(argv[i]) == "--catlist" ||
       string(argv[i]) == "-cl"      || string(argv[i]) == "--cl") {
      if(i+1 < argc) { catlistfile = argv[++i]; }
      else { cerr << "ERROR: -catlist requires a filename argument\n"; return(1); }
    } else if(string(argv[i]) == "-out"    || string(argv[i]) == "--out" ||
              string(argv[i]) == "-output" || string(argv[i]) == "--output" ||
              string(argv[i]) == "-o"      || string(argv[i]) == "--o") {
      if(i+1 < argc) { outfile = argv[++i]; }
      else { cerr << "ERROR: -out requires a filename argument\n"; return(1); }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "--verbose" ||
              string(argv[i]) == "-v"       || string(argv[i]) == "--v") {
      if(i+1 < argc) { verbose = stoi(argv[++i]); }
      else { cerr << "ERROR: -verbose requires an integer argument\n"; return(1); }
    } else {
      cerr << "WARNING: unrecognized argument " << argv[i] << " ignored\n";
    }
    i++;
  }

  if(catlistfile.empty()) {
    cerr << "ERROR: no catalog list file specified (-catlist)\n";
    show_usage(); return(1);
  }
  if(outfile.empty()) {
    cerr << "ERROR: no output file specified (-out)\n";
    show_usage(); return(1);
  }

  cout << "Catalog list file: " << catlistfile << "\n";
  cout << "Output file:       " << outfile << "\n";

  // ------------------------------------------------------------------
  // 1. Read the catlist into a vector of (catalog, colformat) pairs.
  //    This is done upfront so we know the total count before launching
  //    the parallel read region.
  // ------------------------------------------------------------------
  vector<pair<string,string>> catlist;
  {
    ifstream instream1(catlistfile);
    if(!instream1) {
      cerr << "ERROR: cannot open catalog list file " << catlistfile << "\n";
      return(1);
    }
    string cf, colfmt;
    while(instream1 >> cf >> colfmt) {
      if(!cf.empty() && !colfmt.empty())
        catlist.push_back({cf, colfmt});
    }
  }

  if(catlist.empty()) {
    cerr << "ERROR: no valid catalog entries found in " << catlistfile << "\n";
    return(1);
  }

  long k = (long)catlist.size();
  cout << "Found " << k << " catalog(s) in list.\n";

  // ------------------------------------------------------------------
  // 2. Read each catalog in parallel (one thread per catalog).
  //    Each thread writes into its own sub-vector; no shared state is
  //    modified during the parallel region.
  // ------------------------------------------------------------------
  vector<vector<hldet>> catvecs(k);
  vector<int> statuses(k, 0);

  #pragma omp parallel for schedule(dynamic)
  for(long ci = 0; ci < k; ci++) {
    int idcol=IDCOL, mjdcol=MJDCOL, racol=RACOL, deccol=DECCOL;
    int magcol=MAGCOL, bandcol=BANDCOL, obscodecol=OBSCODECOL;
    int trail_len_col=-1, trail_PA_col=-1, sigmag_col=-1;
    int sig_across_col=-1, sig_along_col=-1;
    int known_obj_col=-1, det_qual_col=-1;

    int cfret = parse_colformat(catlist[ci].second,
                                idcol, mjdcol, racol, deccol, magcol,
                                bandcol, obscodecol,
                                trail_len_col, trail_PA_col,
                                sigmag_col, sig_across_col, sig_along_col,
                                known_obj_col, det_qual_col);
    if(cfret < 0) {
      cerr << "ERROR: cannot open colformat file " << catlist[ci].second << "\n";
      statuses[ci] = 1;
      continue;
    }

    // Pass verbose=0 inside the parallel region to avoid interleaved output;
    // forcerun=1 so missing optional fields use defaults rather than aborting.
    int st = read_detection_filemt2(catlist[ci].first,
                                    mjdcol, racol, deccol, magcol,
                                    idcol, bandcol, obscodecol,
                                    trail_len_col, trail_PA_col,
                                    sigmag_col, sig_across_col, sig_along_col,
                                    known_obj_col, det_qual_col,
                                    catvecs[ci], 0 /*verbose*/, 1 /*forcerun*/);
    statuses[ci] = st;
  }

  // Check for errors and print per-catalog counts sequentially.
  long total = 0;
  for(long ci = 0; ci < k; ci++) {
    if(statuses[ci] != 0) {
      cerr << "ERROR: failed to read catalog " << catlist[ci].first
           << " (status=" << statuses[ci] << ")\n";
      return(statuses[ci]);
    }
    cout << "  Read " << catvecs[ci].size()
         << " detections from " << catlist[ci].first << "\n";
    if(verbose >= 1)
      cout << "    colformat: " << catlist[ci].second << "\n";
    total += (long)catvecs[ci].size();
  }
  cout << "Total detections from " << k << " catalog(s): " << total << "\n";

  // ------------------------------------------------------------------
  // 3. Sort each sub-catalog independently (parallel over catalogs).
  //    Telescope data is almost always already time-ordered, so the
  //    is_sorted check avoids O(N log N) work for the common case at
  //    the cost of a single O(N) pass.
  // ------------------------------------------------------------------
  #pragma omp parallel for schedule(dynamic)
  for(long ci = 0; ci < k; ci++) {
    if(!is_sorted(catvecs[ci].begin(), catvecs[ci].end(), early_hldet()))
      sort(catvecs[ci].begin(), catvecs[ci].end(), early_hldet());
  }

  // ------------------------------------------------------------------
  // 4. K-way merge using a min-heap.
  //    Complexity: O(N log k) — far better than O(N log N) for a
  //    global re-sort when each sub-catalog is already sorted.
  //    Elements are std::moved out of the sub-vectors to avoid copying
  //    string members (idstring / band / obscode).
  // ------------------------------------------------------------------
  vector<hldet> alldet;
  alldet.reserve(total);

  // Heap element: (MJD, catalog_index, position_within_catalog)
  using PQElem = tuple<double, long, long>;
  priority_queue<PQElem, vector<PQElem>, greater<PQElem>> pq;

  vector<long> pos(k, 0);
  for(long ci = 0; ci < k; ci++) {
    if(!catvecs[ci].empty())
      pq.push(make_tuple(catvecs[ci][0].MJD, ci, 0L));
  }

  while(!pq.empty()) {
    double mjd; long ci, di;
    tie(mjd, ci, di) = pq.top(); pq.pop();
    alldet.push_back(std::move(catvecs[ci][di]));
    long next = di + 1;
    if(next < (long)catvecs[ci].size())
      pq.push(make_tuple(catvecs[ci][next].MJD, ci, next));
  }

  if(!alldet.empty()) {
    cout << "MJD range: " << fixed << setprecision(6)
         << alldet.front().MJD << " to " << alldet.back().MJD << "\n";
  }

  // ------------------------------------------------------------------
  // 5. Write merged output catalog.
  //    FILE* + fprintf is substantially faster than ofstream with
  //    repeated setprecision/fixed state changes.
  // ------------------------------------------------------------------
  FILE *fp = fopen(outfile.c_str(), "w");
  if(!fp) {
    cerr << "ERROR: cannot open output file " << outfile << "\n";
    return(1);
  }

  fprintf(fp, "#ID,MJD,RA,Dec,Mag,Band,ObsCode,"
              "trail_len,trail_PA,sigmag,sig_across,sig_along,"
              "known_obj,det_qual\n");

  for(long j = 0; j < (long)alldet.size(); j++) {
    const hldet &d = alldet[j];
    fprintf(fp, "%s,%.7f,%.7f,%.7f,%.4f,%s,%s,%.2f,%.2f,%.4f,%.3f,%.3f,%ld,%ld\n",
            d.idstring, d.MJD, d.RA, d.Dec, (double)d.mag,
            d.band, d.obscode,
            (double)d.trail_len, (double)d.trail_PA,
            (double)d.sigmag, (double)d.sig_across, (double)d.sig_along,
            d.known_obj, d.det_qual);
  }
  fclose(fp);

  cout << "Wrote " << alldet.size() << " detections to " << outfile << "\n";
  cout << "Output colformat for use with make_tracklets:\n";
  cout << "  IDCOL 1\n  MJDCOL 2\n  RACOL 3\n  DECCOL 4\n";
  cout << "  MAGCOL 5\n  BANDCOL 6\n  OBSCODECOL 7\n";
  cout << "  TRAILLENCOL 8\n  TRAILPACOL 9\n  SIGMAGCOL 10\n";
  cout << "  SIGACROSSCOL 11\n  SIGALONGCOL 12\n";
  cout << "  KNOWNOBJCOL 13\n  DETQUALCOL 14\n";

  return(0);
}
