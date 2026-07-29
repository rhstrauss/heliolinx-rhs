// label_hldet_mpc.cpp
// Cross-matches a heliolinx detection catalog against MPC 80-column
// observation files to label which catalog detections already exist
// in the MPC database. Uses a 4D k-d tree (time + sky position) for
// fast spatial-temporal matching.
//
// Unlike label_hldet, which requires a pre-computed ephemeris catalog,
// this tool directly reads MPC 80-column observation files and matches
// purely on RA, Dec, and MJD proximity.
//
// Usage:
//   label_hldet_mpc -pairdets pairdets.csv -colformat colformat.txt \
//     -mpcobs T05.txt T08.txt M22.txt W68.txt \
//     -matchrad 2.0 -timerad 5.0 -outfile labeled.csv
//
// The -matchrad is in arcseconds (spatial match tolerance).
// The -timerad is in seconds (temporal match tolerance).
// Both must be satisfied for a match.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define DAY_TO_DEG_CONV 24.0
#define COLS_TO_READ 14

static void show_usage()
{
  cerr << "Usage: label_hldet_mpc -pairdets pairdets_file -colformat column_format_file\n";
  cerr << "       -mpcobs mpc80_file1 [mpc80_file2 ...]\n";
  cerr << "       -matchrad match_arcsec -timerad time_seconds\n";
  cerr << "       -matchobscode 0/1 -outfile outfile\n";
  cerr << "\n";
  cerr << "Required arguments:\n";
  cerr << "  -pairdets   Heliolinx paired detection catalog (CSV)\n";
  cerr << "  -mpcobs     One or more MPC 80-column observation files\n";
  cerr << "  -outfile    Output labeled detection catalog\n";
  cerr << "\n";
  cerr << "Optional arguments:\n";
  cerr << "  -colformat  Column format file for pairdets (default: standard hldet columns)\n";
  cerr << "  -matchrad   Spatial match radius in arcseconds (default: 2.0)\n";
  cerr << "  -timerad    Temporal match tolerance in seconds (default: 5.0)\n";
  cerr << "  -matchobscode  Require observatory code match: 0=no, 1=yes (default: 1)\n";
}

int main(int argc, char *argv[])
{
  ofstream outstream1;
  ifstream instream1;
  string colformatfile, stest;
  vector<hldet> cat_dets = {};
  vector<hldet> mpc_dets = {};
  string pairdets_file, outfile;
  vector<string> mpcobs_files = {};
  double matchrad = 2.0;   // arcseconds
  double timerad = 5.0;    // seconds
  int matchobscode = 1;
  long catnum, mpcnum, catct, mpcct, i;
  catnum = mpcnum = catct = mpcct = i = 0;
  int verbose = 0;
  long kdroot = 0;
  long splitpoint = 0;
  long index = 0;
  long colreadct;
  long nearest = 0;
  double dist = 0.0;
  int colformatfile_set = 0;
  point4d_index onepoint = point4d_index(0, 0, 0, 0, 0);
  point4d_index querypoint = point4d_index(0, 0, 0, 0, 0);
  vector<point4d_index> poolvec;
  KD_point4d_index kdpoint = KD_point4d_index(onepoint, -1, -1, 1, 0);
  vector<KD_point4d_index> kdvec;
  vector<long> indexvec;
  double mjdref = 0.0;
  double timescale = DAY_TO_DEG_CONV;
  int status = 0;

  // Column indices (defaults for standard hldet format)
  int mjdcol = 1, racol = 2, deccol = 3, magcol = 4;
  int idcol = 11, bandcol = 12, obscodecol = 13;
  int trail_len_col = 5, trail_PA_col = 6, sigmag_col = 7;
  int sig_across_col = 8, sig_along_col = 9;
  int known_obj_col = 14, det_qual_col = 15;

  if(argc < 5) {
    show_usage();
    return(1);
  }

  // Parse command-line arguments
  i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-pairdets" || string(argv[i]) == "-pd" || string(argv[i]) == "--pairdets") {
      if(i + 1 < argc) {
        pairdets_file = argv[++i];
        i++;
      } else {
        cerr << "ERROR: -pairdets requires an argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-cf" || string(argv[i]) == "--colformat") {
      if(i + 1 < argc) {
        colformatfile = argv[++i];
        colformatfile_set = 1;
        i++;
      } else {
        cerr << "ERROR: -colformat requires an argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-mpcobs" || string(argv[i]) == "-mpc" || string(argv[i]) == "--mpcobs") {
      // Read all following arguments until next flag or end
      i++;
      while(i < argc && argv[i][0] != '-') {
        mpcobs_files.push_back(string(argv[i]));
        i++;
      }
      if(mpcobs_files.size() == 0) {
        cerr << "ERROR: -mpcobs requires at least one file argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-matchrad" || string(argv[i]) == "-mrad" || string(argv[i]) == "--matchrad") {
      if(i + 1 < argc) {
        matchrad = stod(argv[++i]);
        i++;
      } else {
        cerr << "ERROR: -matchrad requires an argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-timerad" || string(argv[i]) == "-trad" || string(argv[i]) == "--timerad") {
      if(i + 1 < argc) {
        timerad = stod(argv[++i]);
        i++;
      } else {
        cerr << "ERROR: -timerad requires an argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-matchobscode" || string(argv[i]) == "-mobscode" || string(argv[i]) == "--matchobscode") {
      if(i + 1 < argc) {
        matchobscode = stoi(argv[++i]);
        i++;
      } else {
        cerr << "ERROR: -matchobscode requires an argument\n";
        return(1);
      }
    } else if(string(argv[i]) == "-outfile" || string(argv[i]) == "-out" || string(argv[i]) == "--outfile") {
      if(i + 1 < argc) {
        outfile = argv[++i];
        i++;
      } else {
        cerr << "ERROR: -outfile requires an argument\n";
        return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] << "\n";
      i++;
    }
  }

  // Validate required arguments
  if(pairdets_file.empty()) {
    cerr << "ERROR: no pairdets file specified\n";
    show_usage();
    return(1);
  }
  if(mpcobs_files.empty()) {
    cerr << "ERROR: no MPC observation files specified\n";
    show_usage();
    return(1);
  }
  if(outfile.empty()) {
    cerr << "ERROR: no output file specified\n";
    show_usage();
    return(1);
  }

  cout.precision(17);
  cout << "Input pairdets file: " << pairdets_file << "\n";
  cout << "MPC observation files: " << mpcobs_files.size() << " file(s)\n";
  for(i = 0; i < long(mpcobs_files.size()); i++) {
    cout << "  " << mpcobs_files[i] << "\n";
  }
  cout << "Spatial match radius: " << matchrad << " arcsec\n";
  cout << "Temporal match tolerance: " << timerad << " seconds\n";
  cout << "Require observatory code match: " << (matchobscode ? "yes" : "no") << "\n";
  cout << "Output file: " << outfile << "\n";

  // Read column formatting file if supplied
  if(colformatfile_set) {
    instream1.open(colformatfile);
    if(!instream1) {
      cerr << "ERROR: unable to open column format file " << colformatfile << "\n";
      return(1);
    }
    colreadct = 0;
    while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct < COLS_TO_READ) {
      instream1 >> stest;
      if(stest == "MJDCOL") {
        instream1 >> mjdcol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "RACOL") {
        instream1 >> racol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "DECCOL") {
        instream1 >> deccol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "MAGCOL") {
        instream1 >> magcol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "TRAILLENCOL") {
        instream1 >> trail_len_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "TRAILPACOL") {
        instream1 >> trail_PA_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "SIGMAGCOL") {
        instream1 >> sigmag_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "SIGACROSSCOL") {
        instream1 >> sig_across_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "SIGALONGCOL") {
        instream1 >> sig_along_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "IDCOL") {
        instream1 >> idcol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "BANDCOL") {
        instream1 >> bandcol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "OBSCODECOL") {
        instream1 >> obscodecol;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "KNOWNOBJCOL") {
        instream1 >> known_obj_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else if(stest == "DETQUALCOL") {
        instream1 >> det_qual_col;
        if(!instream1.eof() && !instream1.fail() && !instream1.bad()) colreadct++;
      } else {
        cout << "WARNING: unrecognized string " << stest << " in column format file\n";
      }
    }
    instream1.close();
    if(colreadct < COLS_TO_READ) {
      cout << "WARNING: only " << colreadct << " of " << COLS_TO_READ << " expected column specs read from " << colformatfile << "\n";
    }
  }

  // Read the pairdets catalog
  cat_dets = {};
  status = read_detection_filemt2(pairdets_file, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol, trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col, cat_dets, verbose, 1);
  if(status != 0) {
    cerr << "ERROR reading pairdets file " << pairdets_file << " (status " << status << ")\n";
    return(1);
  }
  catnum = cat_dets.size();
  cout << "Read " << catnum << " detections from catalog file " << pairdets_file << "\n";

  // Read all MPC observation files
  mpc_dets = {};
  for(i = 0; i < long(mpcobs_files.size()); i++) {
    long prev_size = mpc_dets.size();
    status = read_detection_file_MPC80(mpcobs_files[i], mpc_dets);
    if(status != 0) {
      cerr << "WARNING: problem reading MPC file " << mpcobs_files[i] << " (status " << status << ")\n";
      cerr << "Continuing with " << mpc_dets.size() << " observations read so far.\n";
    } else {
      cout << "Read " << (mpc_dets.size() - prev_size) << " observations from " << mpcobs_files[i] << "\n";
    }
  }
  mpcnum = mpc_dets.size();
  cout << "Total MPC observations loaded: " << mpcnum << "\n";

  if(mpcnum == 0) {
    cerr << "ERROR: no MPC observations loaded\n";
    return(1);
  }

  // Find the median MJD of catalog detections for the time reference
  {
    vector<double> mjdvec;
    for(catct = 0; catct < catnum; catct++) mjdvec.push_back(cat_dets[catct].MJD);
    mjdref = dmedian(mjdvec);
  }
  cout << "Median MJD of catalog detections: " << mjdref << "\n";

  // Build k-d tree from catalog detections (the larger set).
  // The time dimension uses the same conversion as label_hldet:
  // 1 day = timescale degrees, then projected to radians.
  poolvec = {};
  for(catct = 0; catct < catnum; catct++) {
    onepoint = point4d_index(
      (cat_dets[catct].MJD - mjdref) * timescale / DEGPRAD,
      cos(cat_dets[catct].RA / DEGPRAD) * cos(cat_dets[catct].Dec / DEGPRAD),
      sin(cat_dets[catct].RA / DEGPRAD) * cos(cat_dets[catct].Dec / DEGPRAD),
      sin(cat_dets[catct].Dec / DEGPRAD),
      catct);
    poolvec.push_back(onepoint);
  }
  cout << "Loaded " << poolvec.size() << " catalog points into pool\n";

  // Build the k-d tree
  kdvec = {};
  kdroot = splitpoint = 0;
  splitpoint = medind_4d_index(poolvec, 1);
  kdpoint = KD_point4d_index(poolvec[splitpoint], -1, -1, 1, 0);
  kdvec.push_back(kdpoint);
  kdtree_4d_index(poolvec, 1, splitpoint, kdroot, kdvec);
  cout << "k-d tree constructed with " << kdvec.size() << " nodes\n";

  // Compute the effective search radius for the k-d tree.
  // The k-d tree uses a combined 4D space where time is scaled.
  // We need a search radius large enough to cover both the spatial
  // match radius AND the time tolerance. The effective 4D radius is:
  //   r_4d = sqrt(r_sky^2 + r_time^2)
  // where r_sky is matchrad in radians, and r_time is timerad
  // converted to the same radians-equivalent using the timescale.
  double matchrad_rad = matchrad / ASECPRAD;  // arcsec -> radians
  double timerad_days = timerad / SOLARDAY;    // seconds -> days
  double timerad_rad = timerad_days * timescale / DEGPRAD;  // days -> radians-equiv
  double search_rad = sqrt(matchrad_rad * matchrad_rad + timerad_rad * timerad_rad);

  cout << "Spatial match radius: " << matchrad << " arcsec = " << matchrad_rad << " rad\n";
  cout << "Temporal match radius: " << timerad << " sec = " << timerad_days << " days = " << timerad_rad << " rad-equiv\n";
  cout << "Combined 4D search radius: " << search_rad << " rad\n";

  // Cross-match: query each MPC observation against the catalog k-d tree
  long total_matches = 0;
  long mpc_matched = 0;
  for(mpcct = 0; mpcct < mpcnum; mpcct++) {
    querypoint = point4d_index(
      (mpc_dets[mpcct].MJD - mjdref) * timescale / DEGPRAD,
      cos(mpc_dets[mpcct].RA / DEGPRAD) * cos(mpc_dets[mpcct].Dec / DEGPRAD),
      sin(mpc_dets[mpcct].RA / DEGPRAD) * cos(mpc_dets[mpcct].Dec / DEGPRAD),
      sin(mpc_dets[mpcct].Dec / DEGPRAD),
      mpcct);

    indexvec = {};
    status = kdrange_4d_index(kdvec, querypoint, search_rad, indexvec);

    // Post-filter: check that BOTH spatial and temporal criteria are met
    int matched_any = 0;
    for(long k = 0; k < long(indexvec.size()); k++) {
      index = kdvec[indexvec[k]].point.index;

      // Check time separation
      double dt_sec = fabs(cat_dets[index].MJD - mpc_dets[mpcct].MJD) * SOLARDAY;
      if(dt_sec > timerad) continue;

      // Check angular separation on sky
      double cos_sep = cos(cat_dets[index].RA / DEGPRAD) * cos(cat_dets[index].Dec / DEGPRAD)
                     * cos(mpc_dets[mpcct].RA / DEGPRAD) * cos(mpc_dets[mpcct].Dec / DEGPRAD)
                     + sin(cat_dets[index].RA / DEGPRAD) * cos(cat_dets[index].Dec / DEGPRAD)
                     * sin(mpc_dets[mpcct].RA / DEGPRAD) * cos(mpc_dets[mpcct].Dec / DEGPRAD)
                     + sin(cat_dets[index].Dec / DEGPRAD)
                     * sin(mpc_dets[mpcct].Dec / DEGPRAD);
      if(cos_sep > 1.0) cos_sep = 1.0;
      if(cos_sep < -1.0) cos_sep = -1.0;
      double sep_arcsec = acos(cos_sep) * ASECPRAD;
      if(sep_arcsec > matchrad) continue;

      // Check observatory code match if required
      if(matchobscode) {
        if(string(cat_dets[index].obscode) != string(mpc_dets[mpcct].obscode)) continue;
      }

      // Match found: label the catalog detection with the MPC designation
      stringncopy01(cat_dets[index].idstring, string(mpc_dets[mpcct].idstring), SHORTSTRINGLEN);
      cat_dets[index].known_obj = 999;
      total_matches++;
      matched_any = 1;
    }
    if(matched_any) mpc_matched++;

    if(mpcct % 100000 == 0 && mpcct > 0) {
      cout << "Processed " << mpcct << " / " << mpcnum << " MPC observations, "
           << total_matches << " catalog matches so far\n";
    }
  }

  cout << "\nCross-match complete:\n";
  cout << "  MPC observations processed: " << mpcnum << "\n";
  cout << "  MPC observations with catalog match: " << mpc_matched << "\n";
  cout << "  Total catalog detections labeled: " << total_matches << "\n";

  // Count unique labeled catalog detections
  long labeled_count = 0;
  for(catct = 0; catct < catnum; catct++) {
    if(cat_dets[catct].known_obj == 999) labeled_count++;
  }
  cout << "  Unique catalog detections labeled: " << labeled_count << " / " << catnum
       << " (" << fixed << setprecision(2) << 100.0 * labeled_count / catnum << "%)\n";

  // Write output
  outstream1.open(outfile);
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(catct = 0; catct < catnum; catct++) {
    outstream1 << fixed << setprecision(7) << cat_dets[catct].MJD << ","
               << cat_dets[catct].RA << "," << cat_dets[catct].Dec << ",";
    outstream1 << fixed << setprecision(4) << cat_dets[catct].mag << ",";
    outstream1 << fixed << setprecision(2) << cat_dets[catct].trail_len << ","
               << cat_dets[catct].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << cat_dets[catct].sigmag << ",";
    outstream1 << fixed << setprecision(3) << cat_dets[catct].sig_across << ","
               << cat_dets[catct].sig_along << ",";
    outstream1 << cat_dets[catct].image << "," << cat_dets[catct].idstring << ","
               << cat_dets[catct].band << ",";
    outstream1 << cat_dets[catct].obscode << "," << cat_dets[catct].known_obj << ",";
    outstream1 << cat_dets[catct].det_qual << "," << cat_dets[catct].index << "\n";
  }
  outstream1.close();
  cout << "Wrote " << catnum << " detections to " << outfile << "\n";

  return(0);
}
