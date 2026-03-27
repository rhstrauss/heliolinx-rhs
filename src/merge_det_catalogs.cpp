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
// Usage:
//   merge_det_catalogs -catlist catlist_file -out output_file
//   [-verbose verbosity_level]

#include "solarsyst_dyn_geo01.h"
#include "cmath"

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

int main(int argc, char *argv[])
{
  string catlistfile, outfile;
  int verbose = 0;
  int status = 0;

  if(argc < 5) {
    show_usage();
    return(1);
  }

  int i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-catlist" || string(argv[i]) == "--catlist" ||
       string(argv[i]) == "-cl"      || string(argv[i]) == "--cl") {
      if(i+1 < argc) {
        catlistfile = argv[++i];
      } else {
        cerr << "ERROR: -catlist requires a filename argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-out"    || string(argv[i]) == "--out" ||
              string(argv[i]) == "-output" || string(argv[i]) == "--output" ||
              string(argv[i]) == "-o"      || string(argv[i]) == "--o") {
      if(i+1 < argc) {
        outfile = argv[++i];
      } else {
        cerr << "ERROR: -out requires a filename argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "--verbose" ||
              string(argv[i]) == "-v"       || string(argv[i]) == "--v") {
      if(i+1 < argc) {
        verbose = stoi(argv[++i]);
      } else {
        cerr << "ERROR: -verbose requires an integer argument\n";
        return(1);
      }
    } else {
      cerr << "WARNING: unrecognized argument " << argv[i] << " ignored\n";
    }
    i++;
  }

  if(catlistfile.empty()) {
    cerr << "ERROR: no catalog list file specified (-catlist)\n";
    show_usage();
    return(1);
  }
  if(outfile.empty()) {
    cerr << "ERROR: no output file specified (-out)\n";
    show_usage();
    return(1);
  }

  cout << "Catalog list file: " << catlistfile << "\n";
  cout << "Output file:       " << outfile << "\n";

  // ------------------------------------------------------------------
  // Read the catalog list and accumulate all detections
  // ------------------------------------------------------------------
  ifstream instream1;
  instream1.open(catlistfile);
  if(!instream1) {
    cerr << "ERROR: cannot open catalog list file " << catlistfile << "\n";
    return(1);
  }

  vector<hldet> alldet;
  string catalogfile, colformatfile;
  long total_catalogs = 0;

  while(instream1 >> catalogfile >> colformatfile) {
    if(catalogfile.empty() || colformatfile.empty()) continue;

    cout << "Reading catalog " << catalogfile
         << " with colformat " << colformatfile << "\n";

    // Set column defaults (same defaults as make_tracklets)
    int idcol      = IDCOL;
    int mjdcol     = MJDCOL;
    int racol      = RACOL;
    int deccol     = DECCOL;
    int magcol     = MAGCOL;
    int bandcol    = BANDCOL;
    int obscodecol = OBSCODECOL;
    int trail_len_col  = -1;
    int trail_PA_col   = -1;
    int sigmag_col     = -1;
    int sig_across_col = -1;
    int sig_along_col  = -1;
    int known_obj_col  = -1;
    int det_qual_col   = -1;

    // Parse the colformat file
    ifstream cfstream;
    cfstream.open(colformatfile);
    if(!cfstream) {
      cerr << "ERROR: cannot open colformat file " << colformatfile << "\n";
      return(1);
    }
    int colreadct = 0;
    string stest;
    while(!cfstream.eof() && !cfstream.fail() && !cfstream.bad()
          && colreadct < COLS_TO_READ) {
      cfstream >> stest;
      if(stest == "MJDCOL")      { cfstream >> mjdcol;       colreadct++; }
      else if(stest == "RACOL")  { cfstream >> racol;        colreadct++; }
      else if(stest == "DECCOL") { cfstream >> deccol;       colreadct++; }
      else if(stest == "MAGCOL") { cfstream >> magcol;       colreadct++; }
      else if(stest == "TRAILLENCOL")  { cfstream >> trail_len_col;  colreadct++; }
      else if(stest == "TRAILPACOL")   { cfstream >> trail_PA_col;   colreadct++; }
      else if(stest == "SIGMAGCOL")    { cfstream >> sigmag_col;     colreadct++; }
      else if(stest == "SIGACROSSCOL") { cfstream >> sig_across_col; colreadct++; }
      else if(stest == "SIGALONGCOL")  { cfstream >> sig_along_col;  colreadct++; }
      else if(stest == "IDCOL")        { cfstream >> idcol;          colreadct++; }
      else if(stest == "BANDCOL")      { cfstream >> bandcol;        colreadct++; }
      else if(stest == "OBSCODECOL")   { cfstream >> obscodecol;     colreadct++; }
      else if(stest == "KNOWNOBJCOL")  { cfstream >> known_obj_col;  colreadct++; }
      else if(stest == "DETQUALCOL")   { cfstream >> det_qual_col;   colreadct++; }
      else {
        if(verbose >= 1)
          cout << "WARNING: unrecognized colformat keyword " << stest << "\n";
      }
    }
    cfstream.close();

    if(verbose >= 1) {
      cout << "  Column assignments: IDCOL=" << idcol
           << " MJDCOL=" << mjdcol
           << " RACOL="  << racol
           << " DECCOL=" << deccol
           << " MAGCOL=" << magcol
           << " BANDCOL=" << bandcol
           << " OBSCODECOL=" << obscodecol << "\n";
    }

    // Read detections from this catalog
    vector<hldet> catvec;
    status = read_detection_filemt2(catalogfile, mjdcol, racol, deccol, magcol,
                                    idcol, bandcol, obscodecol,
                                    trail_len_col, trail_PA_col,
                                    sigmag_col, sig_across_col, sig_along_col,
                                    known_obj_col, det_qual_col,
                                    catvec, verbose, 1 /*forcerun*/);
    if(status != 0) {
      cerr << "ERROR: read_detection_filemt2 returned status " << status
           << " for catalog " << catalogfile << "\n";
      return(status);
    }
    cout << "  Read " << catvec.size() << " detections from " << catalogfile << "\n";

    for(auto &d : catvec) alldet.push_back(d);
    total_catalogs++;
  }
  instream1.close();

  if(total_catalogs == 0) {
    cerr << "ERROR: no valid catalog entries found in " << catlistfile << "\n";
    return(1);
  }

  cout << "Total detections from " << total_catalogs << " catalog(s): "
       << alldet.size() << "\n";

  // ------------------------------------------------------------------
  // Sort by MJD
  // ------------------------------------------------------------------
  sort(alldet.begin(), alldet.end(), early_hldet());
  cout << "Sorted " << alldet.size() << " detections by MJD.\n";

  if(alldet.size() > 0) {
    cout << "MJD range: " << fixed << setprecision(6)
         << alldet.front().MJD << " to " << alldet.back().MJD << "\n";
  }

  // ------------------------------------------------------------------
  // Write merged output catalog
  // ------------------------------------------------------------------
  ofstream outstream1;
  outstream1.open(outfile);
  if(!outstream1) {
    cerr << "ERROR: cannot open output file " << outfile << "\n";
    return(1);
  }

  outstream1 << "#ID,MJD,RA,Dec,Mag,Band,ObsCode,"
             << "trail_len,trail_PA,sigmag,sig_across,sig_along,"
             << "known_obj,det_qual\n";

  for(long j = 0; j < long(alldet.size()); j++) {
    const hldet &d = alldet[j];
    outstream1 << d.idstring << ",";
    outstream1 << fixed << setprecision(7) << d.MJD << ",";
    outstream1 << fixed << setprecision(7) << d.RA  << ",";
    outstream1 << fixed << setprecision(7) << d.Dec << ",";
    outstream1 << fixed << setprecision(4) << d.mag  << ",";
    outstream1 << d.band << ",";
    outstream1 << d.obscode << ",";
    outstream1 << fixed << setprecision(2) << d.trail_len << ",";
    outstream1 << fixed << setprecision(2) << d.trail_PA  << ",";
    outstream1 << fixed << setprecision(4) << d.sigmag     << ",";
    outstream1 << fixed << setprecision(3) << d.sig_across << ",";
    outstream1 << fixed << setprecision(3) << d.sig_along  << ",";
    outstream1 << d.known_obj << ",";
    outstream1 << d.det_qual  << "\n";
  }
  outstream1.close();

  cout << "Wrote " << alldet.size() << " detections to " << outfile << "\n";
  cout << "Output colformat for use with make_tracklets:\n";
  cout << "  IDCOL 1\n  MJDCOL 2\n  RACOL 3\n  DECCOL 4\n";
  cout << "  MAGCOL 5\n  BANDCOL 6\n  OBSCODECOL 7\n";
  cout << "  TRAILLENCOL 8\n  TRAILPACOL 9\n  SIGMAGCOL 10\n";
  cout << "  SIGACROSSCOL 11\n  SIGALONGCOL 12\n";
  cout << "  KNOWNOBJCOL 13\n  DETQUALCOL 14\n";

  return(0);
}
