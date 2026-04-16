// split_tracklets_by_time.cpp: March 2026
// Time-partitions the four output files from make_tracklets (or
// make_trailed_tracklets) into N non-overlapping windows, each
// covering at most X days. This allows running heliolinc separately
// on each time slice without re-running make_tracklets.
//
// Partitioning logic:
//   - Each tracklet is assigned to the window that contains its first
//     image (Img1). Because tracklets are intra-night, the second image
//     (Img2) will normally fall in the same window. If the requested
//     window size is smaller than a typical night (~0.6 days) some
//     tracklets may reference images from adjacent windows; a warning
//     is printed but the run continues, including all referenced images
//     in the window's output image file.
//   - Pairdets are included based solely on which tracklets reference
//     them; "orphan" pairdets (detections in no tracklet) are dropped.
//   - Windows with no tracklets produce no output files.
//
// Output naming:
//   The stem is derived from the -imgs filename by stripping the leading
//   "outim_" and the file extension. For example, -imgs outim_foo_01.txt
//   yields stem "foo_01", and window 1 produces:
//     outim_foo_01_split001.txt
//     pairdets_foo_01_split001.csv
//     tracklets_foo_01_split001.csv
//     trk2det_foo_01_split001.csv
//   Override with -outstem to supply a custom stem.
//   Zero-padding uses enough digits to represent the largest window number
//   (minimum 3 digits).
//
// Usage:
//   split_tracklets_by_time
//     -imgs       input_image_file
//     -pairdets   input_paired_detection_file
//     -tracklets  input_tracklet_file
//     -trk2det    input_trk2det_file
//     -window     window_size_days
//     [-outstem   custom_stem]
//     [-verbose   verbosity]

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <iomanip>
#include <sstream>
#include <set>
#include <algorithm>

static void show_usage()
{
  cerr << "Usage: split_tracklets_by_time -imgs imfile -pairdets pairdetfile\n";
  cerr << "  -tracklets trackletfile -trk2det trk2detfile\n";
  cerr << "  -window window_size_days\n";
  cerr << "  [-outstem custom_stem] [-noon_local utc_offset_hours] [-verbose verbosity]\n";
  cerr << "\nDefault output stem is derived from -imgs by stripping \"outim_\" and extension.\n";
  cerr << "Output files per window: outim_{stem}_split{NNN}.txt,\n";
  cerr << "  pairdets_{stem}_split{NNN}.csv, tracklets_{stem}_split{NNN}.csv,\n";
  cerr << "  trk2det_{stem}_split{NNN}.csv\n";
  cerr << "\n-noon_local: snap window boundaries to local noon for the given UTC offset.\n";
  cerr << "  For Hawaii (HST = UTC-10), use -noon_local -10. This places boundaries at\n";
  cerr << "  22:00 UTC (local noon), which falls in daytime gaps and minimizes\n";
  cerr << "  cross-boundary tracklets. Default: off (boundaries start at min MJD).\n";
}

// Zero-padded window label, 1-indexed, minimum width 3.
static string window_label(long w, int width)
{
  ostringstream oss;
  oss << setw(width) << setfill('0') << (w + 1);
  return oss.str();
}

// Extract stem from image filename:
//   strip any leading directory, strip "outim_" prefix, strip extension.
static string stem_from_imfile(const string &imfile)
{
  // Strip directory
  string base = imfile;
  size_t slash = base.rfind('/');
  if(slash != string::npos) base = base.substr(slash + 1);
  // Strip "outim_" prefix
  if(base.size() > 6 && base.substr(0, 6) == "outim_") base = base.substr(6);
  // Strip extension
  size_t dot = base.rfind('.');
  if(dot != string::npos) base = base.substr(0, dot);
  return base;
}

int main(int argc, char *argv[])
{
  string imfile, pairdetfile, trackletfile, trk2detfile;
  string outstem;
  double window_days = 0.0;
  double noon_local_offset = NAN;
  int verbose = 0;
  int status = 0;
  long i;

  if(argc < 11) {
    show_usage();
    return(1);
  }

  i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-imgs" || string(argv[i]) == "--imgs") {
      if(i+1 < argc) { imfile = argv[++i]; i++; }
      else { cerr << "ERROR: -imgs requires a filename\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-pairdets" || string(argv[i]) == "--pairdets") {
      if(i+1 < argc) { pairdetfile = argv[++i]; i++; }
      else { cerr << "ERROR: -pairdets requires a filename\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "--tracklets") {
      if(i+1 < argc) { trackletfile = argv[++i]; i++; }
      else { cerr << "ERROR: -tracklets requires a filename\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "--trk2det") {
      if(i+1 < argc) { trk2detfile = argv[++i]; i++; }
      else { cerr << "ERROR: -trk2det requires a filename\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-window" || string(argv[i]) == "--window") {
      if(i+1 < argc) { window_days = stod(argv[++i]); i++; }
      else { cerr << "ERROR: -window requires a number\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-outstem" || string(argv[i]) == "--outstem") {
      if(i+1 < argc) { outstem = argv[++i]; i++; }
      else { cerr << "ERROR: -outstem requires a string\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-noon_local" || string(argv[i]) == "--noon_local") {
      if(i+1 < argc) { noon_local_offset = stod(argv[++i]); i++; }
      else { cerr << "ERROR: -noon_local requires a UTC offset in hours\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "--verbose") {
      if(i+1 < argc) { verbose = stoi(argv[++i]); i++; }
      else { cerr << "ERROR: -verbose requires a number\n"; show_usage(); return(1); }
    } else {
      cerr << "Warning: unrecognized argument " << argv[i] << "\n";
      i++;
    }
  }

  // Validate required arguments
  if(imfile.empty())       { cerr << "ERROR: -imgs is required\n";      show_usage(); return(1); }
  if(pairdetfile.empty())  { cerr << "ERROR: -pairdets is required\n";  show_usage(); return(1); }
  if(trackletfile.empty()) { cerr << "ERROR: -tracklets is required\n"; show_usage(); return(1); }
  if(trk2detfile.empty())  { cerr << "ERROR: -trk2det is required\n";   show_usage(); return(1); }
  if(window_days <= 0.0)   { cerr << "ERROR: -window must be > 0\n";    show_usage(); return(1); }

  // Derive stem from image filename if not overridden
  if(outstem.empty()) outstem = stem_from_imfile(imfile);
  cout << "Output stem: " << outstem << "\n";

  // Read input files
  vector<hlimage> images;
  status = read_image_file2(imfile, images);
  if(status != 0) { cerr << "ERROR: could not read image file " << imfile << "\n"; return(1); }
  cout << "Read " << images.size() << " images from " << imfile << "\n";

  vector<hldet> pairdets;
  status = read_pairdet_file(pairdetfile, pairdets, verbose);
  if(status != 0) { cerr << "ERROR: could not read pairdet file " << pairdetfile << "\n"; return(1); }
  cout << "Read " << pairdets.size() << " paired detections from " << pairdetfile << "\n";

  vector<tracklet> tracklets;
  status = read_tracklet_file(trackletfile, tracklets, verbose);
  if(status != 0) { cerr << "ERROR: could not read tracklet file " << trackletfile << "\n"; return(1); }
  cout << "Read " << tracklets.size() << " tracklets from " << trackletfile << "\n";

  vector<longpair> trk2det;
  status = read_longpair_file(trk2detfile, trk2det, verbose);
  if(status != 0) { cerr << "ERROR: could not read trk2det file " << trk2detfile << "\n"; return(1); }
  cout << "Read " << trk2det.size() << " trk2det entries from " << trk2detfile << "\n";

  long nimages    = images.size();
  long npairdets  = pairdets.size();
  long ntracklets = tracklets.size();
  long ntrk2det   = trk2det.size();

  if(nimages == 0 || ntracklets == 0) {
    cerr << "ERROR: input files are empty\n";
    return(1);
  }

  // Find MJD range from image file
  double min_mjd = images[0].MJD;
  double max_mjd = images[0].MJD;
  for(i = 1; i < nimages; i++) {
    if(images[i].MJD < min_mjd) min_mjd = images[i].MJD;
    if(images[i].MJD > max_mjd) max_mjd = images[i].MJD;
  }
  cout << "MJD range: " << fixed << setprecision(5) << min_mjd << " to " << max_mjd << "\n";

  // Compute window origin: either min_mjd or snapped to local noon
  double win_origin = min_mjd;
  if(!isnan(noon_local_offset)) {
    double noon_utc_frac = fmod((12.0 - noon_local_offset) / 24.0, 1.0);
    if(noon_utc_frac < 0.0) noon_utc_frac += 1.0;
    win_origin = floor(min_mjd) + noon_utc_frac;
    if(win_origin > min_mjd) win_origin -= 1.0;
    cout << "Noon-local snapping: UTC offset = " << noon_local_offset
         << " h, noon at MJD fractional " << fixed << setprecision(5) << noon_utc_frac
         << ", window origin = " << win_origin << "\n";
  }

  // Compute number of windows
  long numwin = (long)ceil((max_mjd - win_origin) / window_days);
  if(numwin < 1) numwin = 1;
  cout << "Splitting into up to " << numwin << " windows of " << window_days << " days each\n";

  // Assign each image to a window
  vector<long> img_win(nimages);
  for(i = 0; i < nimages; i++) {
    long w = (long)floor((images[i].MJD - win_origin) / window_days);
    if(w < 0) w = 0;
    if(w >= numwin) w = numwin - 1;
    img_win[i] = w;
  }

  // Assign each tracklet to a window based on its first image (Img1)
  vector<long> trk_win(ntracklets);
  long cross_boundary_count = 0;
  for(long t = 0; t < ntracklets; t++) {
    long img1 = tracklets[t].Img1;
    if(img1 < 0 || img1 >= nimages) {
      cerr << "ERROR: tracklet " << t << " has out-of-range Img1=" << img1 << "\n";
      return(1);
    }
    trk_win[t] = img_win[img1];
    long img2 = tracklets[t].Img2;
    if(img2 < 0 || img2 >= nimages) {
      cerr << "ERROR: tracklet " << t << " has out-of-range Img2=" << img2 << "\n";
      return(1);
    }
    if(img_win[img2] != trk_win[t]) {
      cross_boundary_count++;
      if(verbose >= 1) {
        cerr << "Warning: tracklet " << t << " spans window boundary (Img1 in window "
             << trk_win[t] << ", Img2 in window " << img_win[img2] << ")\n";
      }
    }
  }
  if(cross_boundary_count > 0) {
    cerr << "Warning: " << cross_boundary_count << " of " << ntracklets
         << " tracklets span a window boundary. Consider -noon_local or a larger -window.\n";
  }

  // Partition trk2det entries by window (single pass)
  vector<vector<longpair>> win_trk2det(numwin);
  for(i = 0; i < ntrk2det; i++) {
    long t = trk2det[i].i1;
    long d = trk2det[i].i2;
    if(t < 0 || t >= ntracklets) {
      cerr << "ERROR: trk2det entry " << i << " has out-of-range tracklet index " << t << "\n";
      return(1);
    }
    if(d < 0 || d >= npairdets) {
      cerr << "ERROR: trk2det entry " << i << " has out-of-range pairdet index " << d << "\n";
      return(1);
    }
    win_trk2det[trk_win[t]].push_back(trk2det[i]);
  }

  // Zero-pad width for window labels
  int label_width = 3;
  {
    long tmp = numwin;
    int w = 0;
    while(tmp > 0) { tmp /= 10; w++; }
    if(w > label_width) label_width = w;
  }

  long windows_written = 0;

  for(long w = 0; w < numwin; w++) {
    // Collect tracklets in this window
    vector<long> win_trks;
    for(long t = 0; t < ntracklets; t++) {
      if(trk_win[t] == w) win_trks.push_back(t);
    }

    if(win_trks.empty()) {
      if(verbose >= 1) cout << "Window " << (w+1) << ": no tracklets, skipping\n";
      continue;
    }

    // Collect pairdet old-indices needed by this window's tracklets
    set<long> det_set;
    for(auto& pair : win_trk2det[w]) det_set.insert(pair.i2);
    vector<long> win_dets(det_set.begin(), det_set.end()); // sorted

    // Collect image old-indices needed: from pairdets and tracklet Img1/Img2
    set<long> img_set;
    for(long d : win_dets) {
      long img = pairdets[d].image;
      if(img < 0 || img >= nimages) {
        cerr << "ERROR: pairdet " << d << " has out-of-range image index " << img << "\n";
        return(1);
      }
      img_set.insert(img);
    }
    for(long t : win_trks) {
      img_set.insert(tracklets[t].Img1);
      img_set.insert(tracklets[t].Img2);
    }
    vector<long> win_imgs(img_set.begin(), img_set.end()); // sorted

    // Build old -> new local index maps
    vector<long> img_old2new(nimages, -1);
    for(long ni = 0; ni < (long)win_imgs.size(); ni++) img_old2new[win_imgs[ni]] = ni;

    vector<long> det_old2new(npairdets, -1);
    for(long nd = 0; nd < (long)win_dets.size(); nd++) det_old2new[win_dets[nd]] = nd;

    vector<long> trk_old2new(ntracklets, -1);
    for(long nt = 0; nt < (long)win_trks.size(); nt++) trk_old2new[win_trks[nt]] = nt;

    string label = window_label(w, label_width);
    string suffix = "_split" + label;

    // --- Write image file ---
    string out_imfile = "outim_" + outstem + suffix + ".txt";
    ofstream out1(out_imfile);
    if(!out1) { cerr << "ERROR: could not open " << out_imfile << " for writing\n"; return(1); }
    for(long old_i : win_imgs) {
      const hlimage& im = images[old_i];
      out1 << fixed << setprecision(8) << im.MJD << " " << im.RA << " " << im.Dec
           << " " << im.obscode << " "
           << fixed << setprecision(1) << im.X << " " << im.Y << " " << im.Z << " "
           << fixed << setprecision(4) << im.VX << " " << im.VY << " " << im.VZ << " "
           << im.startind << " " << im.endind << " " << im.exptime << "\n";
    }
    out1.close();

    // --- Write pairdets file ---
    string out_pairdetfile = "pairdets_" + outstem + suffix + ".csv";
    ofstream out2(out_pairdetfile);
    if(!out2) { cerr << "ERROR: could not open " << out_pairdetfile << " for writing\n"; return(1); }
    out2 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
    for(long old_d : win_dets) {
      const hldet& d = pairdets[old_d];
      long new_img = img_old2new[d.image];
      if(new_img < 0) {
        cerr << "ERROR: pairdet " << old_d << " image index " << d.image
             << " not in window " << w << " image set\n";
        return(1);
      }
      out2 << fixed << setprecision(7) << d.MJD << "," << d.RA << "," << d.Dec << ","
           << fixed << setprecision(4) << d.mag << ","
           << fixed << setprecision(2) << d.trail_len << "," << d.trail_PA << ","
           << fixed << setprecision(4) << d.sigmag << ","
           << fixed << setprecision(3) << d.sig_across << "," << d.sig_along << ","
           << new_img << "," << d.idstring << "," << d.band << ","
           << d.obscode << "," << d.known_obj << ","
           << d.det_qual << "," << d.index << "\n";
    }
    out2.close();

    // --- Write tracklets file ---
    string out_trkfile = "tracklets_" + outstem + suffix + ".csv";
    ofstream out3(out_trkfile);
    if(!out3) { cerr << "ERROR: could not open " << out_trkfile << " for writing\n"; return(1); }
    out3 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
    for(long nt = 0; nt < (long)win_trks.size(); nt++) {
      long old_t = win_trks[nt];
      const tracklet& trk = tracklets[old_t];
      long new_img1 = img_old2new[trk.Img1];
      long new_img2 = img_old2new[trk.Img2];
      if(new_img1 < 0 || new_img2 < 0) {
        cerr << "ERROR: tracklet " << old_t << " image indices not in window " << w << "\n";
        return(1);
      }
      out3 << fixed << setprecision(7)
           << new_img1 << "," << trk.RA1 << "," << trk.Dec1 << ","
           << new_img2 << "," << trk.RA2 << "," << trk.Dec2 << ","
           << trk.npts << "," << nt << "\n";
    }
    out3.close();

    // --- Write trk2det file ---
    string out_trk2detfile = "trk2det_" + outstem + suffix + ".csv";
    ofstream out4(out_trk2detfile);
    if(!out4) { cerr << "ERROR: could not open " << out_trk2detfile << " for writing\n"; return(1); }
    out4 << "#trk_ID,detnum\n";
    for(auto& pair : win_trk2det[w]) {
      long new_t = trk_old2new[pair.i1];
      long new_d = det_old2new[pair.i2];
      if(new_t < 0 || new_d < 0) {
        cerr << "ERROR: trk2det entry (" << pair.i1 << "," << pair.i2
             << ") remapping failed in window " << w << "\n";
        return(1);
      }
      out4 << new_t << "," << new_d << "\n";
    }
    out4.close();

    windows_written++;
    cout << "Window " << (w+1) << " [" << fixed << setprecision(3)
         << (win_origin + w * window_days) << " to "
         << (win_origin + (w+1) * window_days) << "]: "
         << win_imgs.size() << " images, "
         << win_dets.size() << " pairdets, "
         << win_trks.size() << " tracklets, "
         << win_trk2det[w].size() << " trk2det entries -> "
         << label << "\n";
  }

  cout << "Wrote " << windows_written << " time windows.\n";
  return(0);
}
