// split_tracklets_by_time.cpp: March 2026
// Time-partitions the four output files from make_tracklets (or
// make_trailed_tracklets) into N windows, each covering at most X days,
// with optional overlap between consecutive windows.
//
// Partitioning logic:
//   - Window w covers the half-open interval
//       [win_origin + w*stride,  win_origin + w*stride + window_days)
//     where stride = window_days - overlap_days.
//   - When overlap_days == 0 (the default) this reproduces the original
//     non-overlapping behaviour exactly.
//   - When overlap_days > 0, consecutive windows share time.  A tracklet
//     whose Img1 MJD falls inside multiple windows is emitted into every
//     one of them with locally renumbered indices; heliolinc therefore
//     processes it in each overlapping run.
//   - Each tracklet is assigned to every window whose time range contains
//     its Img1 MJD.  Because tracklets are intra-night, Img2 will
//     normally also fall inside the same window(s).  If a tracklet's Img2
//     MJD lies outside a window that owns its Img1 (possible when
//     stride < one night or when -noon_local is not used), a warning is
//     printed but the run continues, including all referenced images in
//     that window's output image file.
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
//     [-overlap   overlap_days]          (default 0; must be >= 0 and < window)
//     [-outstem   custom_stem]
//     [-noon_local utc_offset_hours]
//     [-verbose   verbosity]

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <algorithm>

static void show_usage()
{
  cerr << "Usage: split_tracklets_by_time -imgs imfile -pairdets pairdetfile\n";
  cerr << "  -tracklets trackletfile -trk2det trk2detfile\n";
  cerr << "  -window window_size_days\n";
  cerr << "  [-overlap overlap_days]  (default 0.0; must be >= 0 and < window)\n";
  cerr << "  [-outstem custom_stem] [-noon_local utc_offset_hours] [-verbose verbosity]\n";
  cerr << "\nDefault output stem is derived from -imgs by stripping \"outim_\" and extension.\n";
  cerr << "Output files per window: outim_{stem}_split{NNN}.txt,\n";
  cerr << "  pairdets_{stem}_split{NNN}.csv, tracklets_{stem}_split{NNN}.csv,\n";
  cerr << "  trk2det_{stem}_split{NNN}.csv\n";
  cerr << "\n-overlap: days shared between consecutive windows. With -window 7 -overlap 5,\n";
  cerr << "  the stride is 2 days, so windows start at t0, t0+2, t0+4, ... and each\n";
  cerr << "  covers 7 days. A tracklet that falls in multiple windows is emitted into all.\n";
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

// 8 MB write buffer per stream — materially faster for multi-GB CSV output.
static const int WBUF_SIZE = 8 * 1024 * 1024;

int main(int argc, char *argv[])
{
  // Detach C stdio from C++ streams; we never mix them.
  ios::sync_with_stdio(false);

  string imfile, pairdetfile, trackletfile, trk2detfile;
  string outstem;
  double window_days  = 0.0;
  double overlap_days = 0.0;
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
    } else if(string(argv[i]) == "-overlap" || string(argv[i]) == "--overlap") {
      if(i+1 < argc) { overlap_days = stod(argv[++i]); i++; }
      else { cerr << "ERROR: -overlap requires a number\n"; show_usage(); return(1); }
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
  if(overlap_days < 0.0)   { cerr << "ERROR: -overlap must be >= 0\n";  show_usage(); return(1); }
  if(overlap_days >= window_days) {
    cerr << "ERROR: -overlap (" << overlap_days
         << ") must be strictly less than -window (" << window_days << ")\n";
    show_usage(); return(1);
  }

  double stride = window_days - overlap_days;

  // Derive stem from image filename if not overridden
  if(outstem.empty()) outstem = stem_from_imfile(imfile);
  cout << "Output stem: " << outstem << "\n";
  if(overlap_days > 0.0) {
    cout << "Overlap mode: window=" << window_days << " d, overlap=" << overlap_days
         << " d, stride=" << stride << " d\n";
  }

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

  // Compute number of windows.
  // With stride < window, the last window must still reach max_mjd.
  // Formula: window w covers [win_origin + w*stride, win_origin + w*stride + window_days).
  // We need win_origin + w*stride <= max_mjd, solved for w:
  //   w_max = ceil((max_mjd - win_origin - window_days) / stride) + 1
  // (clamp to 1 for the degenerate case).
  long numwin;
  if(stride >= (max_mjd - win_origin)) {
    numwin = 1;
  } else {
    numwin = (long)ceil((max_mjd - win_origin - window_days) / stride) + 1;
    if(numwin < 1) numwin = 1;
  }
  cout << "Splitting into up to " << numwin << " windows of " << window_days
       << " days each (stride " << stride << " d)\n";

  // Build per-image MJD lookup (images vector is already loaded).
  // img_mjd[i] == images[i].MJD — we read it inline below for clarity.

  // Assign each tracklet to its windows (single pass, O(ntracklets * avg_windows_per_trk)).
  // For the no-overlap case each tracklet hits exactly one window; for overlap it may hit 2.
  //
  // win_trks[w] holds the global (input) tracklet indices assigned to window w.
  vector<vector<long>> win_trks(numwin);
  // Pre-size: in the non-overlap case every tracklet lands in exactly one window; reserve
  // a proportional share so we don't thrash the allocator.
  {
    long avg_per_win = (ntracklets + numwin - 1) / numwin;
    for(long w = 0; w < numwin; w++) win_trks[w].reserve(avg_per_win);
  }

  long cross_boundary_count = 0;

  for(long t = 0; t < ntracklets; t++) {
    long img1 = tracklets[t].Img1;
    if(img1 < 0 || img1 >= nimages) {
      cerr << "ERROR: tracklet " << t << " has out-of-range Img1=" << img1 << "\n";
      return(1);
    }
    long img2 = tracklets[t].Img2;
    if(img2 < 0 || img2 >= nimages) {
      cerr << "ERROR: tracklet " << t << " has out-of-range Img2=" << img2 << "\n";
      return(1);
    }
    double mjd1 = images[img1].MJD;
    double mjd2 = images[img2].MJD;

    // Find first window whose start <= mjd1 < start + window_days.
    // Window w starts at win_origin + w * stride.
    // mjd1 >= win_origin + w*stride  =>  w <= (mjd1 - win_origin) / stride
    // mjd1 <  win_origin + w*stride + window_days  =>  w > (mjd1 - win_origin - window_days) / stride
    double offset = mjd1 - win_origin;
    long w_lo = (offset >= window_days)
                  ? (long)floor((offset - window_days) / stride) + 1
                  : 0;
    long w_hi = (offset >= 0.0)
                  ? (long)floor(offset / stride)
                  : -1;
    if(w_lo < 0)    w_lo = 0;
    if(w_hi >= numwin) w_hi = numwin - 1;

    if(w_lo > w_hi) {
      // Img1 is before win_origin or after all windows — clamp to nearest.
      // This mirrors the original clamping behaviour.
      if(offset < 0.0) {
        win_trks[0].push_back(t);
      } else {
        win_trks[numwin - 1].push_back(t);
      }
    } else {
      for(long w = w_lo; w <= w_hi; w++) {
        win_trks[w].push_back(t);
        // Cross-boundary check: does Img2 fall outside this particular window?
        double wstart = win_origin + w * stride;
        double wend   = wstart + window_days;
        if(mjd2 < wstart || mjd2 >= wend) {
          cross_boundary_count++;
          if(verbose >= 1) {
            cerr << "Warning: tracklet " << t
                 << " spans window boundary in window " << (w+1)
                 << " (Img1 MJD=" << fixed << setprecision(5) << mjd1
                 << " inside [" << wstart << "," << wend
                 << "), Img2 MJD=" << mjd2 << " outside)\n";
          }
        }
      }
    }
  }

  if(cross_boundary_count > 0) {
    cerr << "Warning: " << cross_boundary_count
         << " (window, tracklet) pairs have Img2 outside the window's time range."
         << " Consider -noon_local or a larger -window.\n";
  }

  // Partition trk2det entries by window.
  // A tracklet that appears in multiple windows contributes its trk2det rows to each.
  // Build win_trk2det[w] = trk2det rows for tracklets in win_trks[w].
  //
  // First, build a reverse map: trk_index -> list of windows it belongs to.
  // For the non-overlap case this is trivially one window per tracklet.
  // For overlap it may be two.  We use a flat per-tracklet structure.
  //
  // Simpler approach that is still O(ntrk2det * avg_wins_per_trk):
  // build a lookup trk2wins[t] and iterate trk2det once.
  vector<vector<long>> trk2wins(ntracklets);
  for(long w = 0; w < numwin; w++) {
    for(long t : win_trks[w]) {
      trk2wins[t].push_back(w);
    }
  }

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
    for(long w : trk2wins[t]) {
      win_trk2det[w].push_back(trk2det[i]);
    }
  }
  // Free the reverse map now that we don't need it.
  { vector<vector<long>>().swap(trk2wins); }

  // Zero-pad width for window labels
  int label_width = 3;
  {
    long tmp = numwin;
    int  ww  = 0;
    while(tmp > 0) { tmp /= 10; ww++; }
    if(ww > label_width) label_width = ww;
  }

  long windows_written = 0;

  // Per-window write buffers (one set of four char arrays, reused each iteration).
  vector<char> wbuf1(WBUF_SIZE), wbuf2(WBUF_SIZE), wbuf3(WBUF_SIZE), wbuf4(WBUF_SIZE);

  for(long w = 0; w < numwin; w++) {
    const vector<long>& wtrks = win_trks[w];

    if(wtrks.empty()) {
      if(verbose >= 1) cout << "Window " << (w+1) << ": no tracklets, skipping\n";
      continue;
    }

    double wstart = win_origin + w * stride;
    double wend   = wstart + window_days;

    // Collect pairdet old-indices needed by this window's tracklets.
    // Use an unordered_map for old->new mapping to avoid dense O(N) allocations
    // (for the PS1 scale, a dense vector of npairdets longs would be ~3 GB per window).
    // We gather unique det indices from win_trk2det[w] and sort them for reproducible output.
    vector<long> win_dets;
    win_dets.reserve(win_trk2det[w].size()); // upper bound
    for(auto& pr : win_trk2det[w]) win_dets.push_back(pr.i2);
    sort(win_dets.begin(), win_dets.end());
    win_dets.erase(unique(win_dets.begin(), win_dets.end()), win_dets.end());

    // Collect image old-indices from pairdets and tracklet Img1/Img2.
    vector<long> win_imgs;
    win_imgs.reserve(win_dets.size() + 2 * wtrks.size()); // rough upper bound
    for(long old_d : win_dets) {
      long img = pairdets[old_d].image;
      if(img < 0 || img >= nimages) {
        cerr << "ERROR: pairdet " << old_d << " has out-of-range image index " << img << "\n";
        return(1);
      }
      win_imgs.push_back(img);
    }
    for(long old_t : wtrks) {
      win_imgs.push_back(tracklets[old_t].Img1);
      win_imgs.push_back(tracklets[old_t].Img2);
    }
    sort(win_imgs.begin(), win_imgs.end());
    win_imgs.erase(unique(win_imgs.begin(), win_imgs.end()), win_imgs.end());

    // Build old->new maps using unordered_map.
    unordered_map<long,long> img_old2new, det_old2new, trk_old2new;
    img_old2new.reserve(win_imgs.size() * 2);
    det_old2new.reserve(win_dets.size() * 2);
    trk_old2new.reserve(wtrks.size() * 2);

    for(long ni = 0; ni < (long)win_imgs.size(); ni++) img_old2new[win_imgs[ni]] = ni;
    for(long nd = 0; nd < (long)win_dets.size(); nd++) det_old2new[win_dets[nd]] = nd;
    for(long nt = 0; nt < (long)wtrks.size();   nt++) trk_old2new[wtrks[nt]]    = nt;

    string label  = window_label(w, label_width);
    string suffix = "_split" + label;

    // --- Write image file ---
    string out_imfile = "outim_" + outstem + suffix + ".txt";
    ofstream out1(out_imfile);
    if(!out1) { cerr << "ERROR: could not open " << out_imfile << " for writing\n"; return(1); }
    out1.rdbuf()->pubsetbuf(wbuf1.data(), WBUF_SIZE);
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
    out2.rdbuf()->pubsetbuf(wbuf2.data(), WBUF_SIZE);
    out2 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
    for(long old_d : win_dets) {
      const hldet& d = pairdets[old_d];
      long new_img = img_old2new.at(d.image);
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
    out3.rdbuf()->pubsetbuf(wbuf3.data(), WBUF_SIZE);
    out3 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
    for(long nt = 0; nt < (long)wtrks.size(); nt++) {
      long old_t = wtrks[nt];
      const tracklet& trk = tracklets[old_t];
      long new_img1 = img_old2new.at(trk.Img1);
      long new_img2 = img_old2new.at(trk.Img2);
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
    out4.rdbuf()->pubsetbuf(wbuf4.data(), WBUF_SIZE);
    out4 << "#trk_ID,detnum\n";
    for(auto& pr : win_trk2det[w]) {
      long new_t = trk_old2new.at(pr.i1);
      long new_d = det_old2new.at(pr.i2);
      out4 << new_t << "," << new_d << "\n";
    }
    out4.close();

    windows_written++;
    cout << "Window " << (w+1) << " [" << fixed << setprecision(3)
         << wstart << " to " << wend << "]: "
         << win_imgs.size() << " images, "
         << win_dets.size() << " pairdets, "
         << wtrks.size() << " tracklets, "
         << win_trk2det[w].size() << " trk2det entries -> "
         << label << "\n";
  }

  cout << "Wrote " << windows_written << " time windows.\n";
  return(0);
}
