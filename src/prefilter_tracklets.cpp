// prefilter_tracklets.cpp: catalog-quality pre-filter for helio_highgrade input files.
// Loads pairdets/tracklets/trk2det/imfile and applies configurable quality filters,
// writing filtered output files suitable for helio_highgrade or heliolinc.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

static void show_usage()
{
  cerr << "Usage: prefilter_tracklets -imgs imfile -pairdets pairdetfile -tracklets trackletfile -trk2det trk2detfile\n";
  cerr << "  [-minrate min_ang_rate_deg_per_day]    default 0.0 (disabled)\n";
  cerr << "  [-maxrate max_ang_rate_deg_per_day]    default 0.0 (disabled)\n";
  cerr << "  [-mindetqual min_det_qual]             default 0 (no filter)\n";
  cerr << "  [-maxdmag max_mag_diff_in_tracklet]    default 0.0 (disabled)\n";
  cerr << "  [-maximdensity max_dets_per_image]     default 0 (disabled)\n";
  cerr << "  [-excludeknown]                        flag: exclude detections with known_obj > 0\n";
  cerr << "  [-maxtraillen max_trail_len_arcsec]    default 0.0 (disabled)\n";
  cerr << "  [-minnights min_distinct_nights]       default 0 (disabled)\n";
  cerr << "  [-nightstep nightstep_days]            default NIGHTSTEP\n";
  cerr << "  [-outpairdets outpairdetfile]          default \"prefiltered_pairdets.csv\"\n";
  cerr << "  [-outtracklets outtrkfile]             default \"prefiltered_tracklets.csv\"\n";
  cerr << "  [-outtrk2det outtrk2detfile]           default \"prefiltered_trk2det.csv\"\n";
  cerr << "  [-verbose verbosity]                   default 0\n";
}

int main(int argc, char *argv[])
{
  cout.precision(17);

  // Input/output filenames
  string imfile, pairdetfile, trackletfile, trk2detfile;
  string outpairdetfile = "prefiltered_pairdets.csv";
  string outtrkfile     = "prefiltered_tracklets.csv";
  string outtrk2detfile = "prefiltered_trk2det.csv";

  // Filter parameters
  double minrate        = 0.0;  // deg/day, 0 = disabled
  double maxrate        = 0.0;  // deg/day, 0 = disabled
  int    mindetqual     = 0;    // 0 = no filter
  double maxdmag        = 0.0;  // 0 = disabled
  long   maximdensity   = 0;    // 0 = disabled
  int    excludeknown   = 0;    // flag
  double maxtraillen    = 0.0;  // arcsec, 0 = disabled
  int    minnights      = 0;    // 0 = disabled
  double nightstep      = NIGHTSTEP;
  int    verbose        = 0;

  if(argc < 9) {
    show_usage();
    return(1);
  }

  // Parse arguments
  int i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-imgs") {
      if(i+1 < argc) { imfile = argv[++i]; i++; }
      else { cerr << "-imgs requires an argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-pairdets") {
      if(i+1 < argc) { pairdetfile = argv[++i]; i++; }
      else { cerr << "-pairdets requires an argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-tracklets") {
      if(i+1 < argc) { trackletfile = argv[++i]; i++; }
      else { cerr << "-tracklets requires an argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-trk2det") {
      if(i+1 < argc) { trk2detfile = argv[++i]; i++; }
      else { cerr << "-trk2det requires an argument\n"; show_usage(); return(1); }
    } else if(string(argv[i]) == "-minrate") {
      if(i+1 < argc) { minrate = stod(argv[++i]); i++; }
      else { cerr << "-minrate requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-maxrate") {
      if(i+1 < argc) { maxrate = stod(argv[++i]); i++; }
      else { cerr << "-maxrate requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-mindetqual") {
      if(i+1 < argc) { mindetqual = stoi(argv[++i]); i++; }
      else { cerr << "-mindetqual requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-maxdmag") {
      if(i+1 < argc) { maxdmag = stod(argv[++i]); i++; }
      else { cerr << "-maxdmag requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-maximdensity") {
      if(i+1 < argc) { maximdensity = stol(argv[++i]); i++; }
      else { cerr << "-maximdensity requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-excludeknown") {
      excludeknown = 1; i++;
    } else if(string(argv[i]) == "-maxtraillen") {
      if(i+1 < argc) { maxtraillen = stod(argv[++i]); i++; }
      else { cerr << "-maxtraillen requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-minnights") {
      if(i+1 < argc) { minnights = stoi(argv[++i]); i++; }
      else { cerr << "-minnights requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-nightstep") {
      if(i+1 < argc) { nightstep = stod(argv[++i]); i++; }
      else { cerr << "-nightstep requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-outpairdets") {
      if(i+1 < argc) { outpairdetfile = argv[++i]; i++; }
      else { cerr << "-outpairdets requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-outtracklets") {
      if(i+1 < argc) { outtrkfile = argv[++i]; i++; }
      else { cerr << "-outtracklets requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-outtrk2det") {
      if(i+1 < argc) { outtrk2detfile = argv[++i]; i++; }
      else { cerr << "-outtrk2det requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-verbose") {
      if(i+1 < argc) { verbose = stoi(argv[++i]); i++; }
      else { cerr << "-verbose requires an argument\n"; return(1); }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  if(imfile.empty() || pairdetfile.empty() || trackletfile.empty() || trk2detfile.empty()) {
    cerr << "ERROR: -imgs, -pairdets, -tracklets, and -trk2det are all required\n";
    show_usage();
    return(1);
  }

  // Read input files
  vector <hlimage> image_log;
  int status = read_image_file2(imfile, image_log);
  if(status != 0) {
    cerr << "ERROR: could not read image file " << imfile << " (status=" << status << ")\n";
    return(1);
  }
  cout << "Read " << image_log.size() << " images from " << imfile << "\n";

  vector <hldet> detvec;
  status = read_pairdet_file(pairdetfile, detvec, verbose);
  if(status != 0) {
    cerr << "ERROR: could not read pairdets file " << pairdetfile << " (status=" << status << ")\n";
    return(1);
  }
  cout << "Read " << detvec.size() << " detections from " << pairdetfile << "\n";

  vector <tracklet> tracklets;
  status = read_tracklet_file(trackletfile, tracklets, verbose);
  if(status != 0) {
    cerr << "ERROR: could not read tracklet file " << trackletfile << " (status=" << status << ")\n";
    return(1);
  }
  cout << "Read " << tracklets.size() << " tracklets from " << trackletfile << "\n";

  vector <longpair> trk2det;
  status = read_longpair_file(trk2detfile, trk2det, verbose);
  if(status != 0) {
    cerr << "ERROR: could not read trk2det file " << trk2detfile << " (status=" << status << ")\n";
    return(1);
  }
  cout << "Read " << trk2det.size() << " trk2det entries from " << trk2detfile << "\n";

  long imnum    = (long)image_log.size();
  long detnum   = (long)detvec.size();
  long pairnum  = (long)tracklets.size();

  // ------------------------------------------------------------------
  // Pass 1: Compute per-image detection count (for maximdensity filter)
  // ------------------------------------------------------------------
  vector<long> image_detcount(imnum, 0);
  if(maximdensity > 0) {
    for(long d = 0; d < detnum; d++) {
      long im = detvec[d].image;
      if(im >= 0 && im < imnum) image_detcount[im]++;
    }
    if(verbose >= 1) {
      cout << "Computed per-image detection counts for maximdensity filter (threshold=" << maximdensity << ")\n";
    }
  }

  // ------------------------------------------------------------------
  // minnights filter: count distinct nights in the dataset
  // ------------------------------------------------------------------
  int total_distinct_nights = 0;
  if(minnights > 0) {
    // Use floor(MJD / nightstep) as the night identifier
    set<long> night_set;
    for(long im = 0; im < imnum; im++) {
      long night_id = (long)floor(image_log[im].MJD / nightstep);
      night_set.insert(night_id);
    }
    total_distinct_nights = (int)night_set.size();
    cout << "Total distinct nights in dataset: " << total_distinct_nights << "\n";
    if(total_distinct_nights < minnights) {
      cout << "WARNING: dataset has only " << total_distinct_nights << " distinct nights, but minnights=" << minnights << ". All tracklets will be rejected by minnights filter.\n";
    }
  }

  // ------------------------------------------------------------------
  // Build O(1) offset index into trk2det for fast per-tracklet lookup
  // ------------------------------------------------------------------
  vector<long> t2d_offsets(pairnum + 1, (long)trk2det.size());
  for(long k = (long)trk2det.size() - 1; k >= 0; k--) {
    long trknum = trk2det[k].i1;
    if(trknum >= 0 && trknum < pairnum) t2d_offsets[trknum] = k;
  }

  // ------------------------------------------------------------------
  // Pass 2: Mark tracklets to keep
  // ------------------------------------------------------------------
  vector<bool> keep_tracklet(pairnum, true);

  // Filter reason counters
  long n_rate_low = 0, n_rate_high = 0, n_detqual = 0, n_dmag = 0;
  long n_imdensity = 0, n_known = 0, n_trail = 0, n_nights = 0;

  for(long t = 0; t < pairnum; t++) {
    if(!keep_tracklet[t]) continue;

    long i1 = tracklets[t].Img1;
    long i2 = tracklets[t].Img2;

    // --- Rate filter ---
    if(minrate > 0.0 || maxrate > 0.0) {
      double dra = tracklets[t].RA2 - tracklets[t].RA1;
      if(dra > 180.0)  dra -= 360.0;
      if(dra < -180.0) dra += 360.0;
      double ddec = tracklets[t].Dec2 - tracklets[t].Dec1;
      double meandec_rad = 0.5*(tracklets[t].Dec1 + tracklets[t].Dec2) * M_PI / 180.0;
      double cosdec = cos(meandec_rad);
      double sep_deg = sqrt(dra*cosdec*dra*cosdec + ddec*ddec);
      double dt = image_log[i2].MJD - image_log[i1].MJD;
      double rate = (dt > 0.0) ? sep_deg / dt : 0.0;
      if(minrate > 0.0 && rate < minrate) {
        keep_tracklet[t] = false;
        n_rate_low++;
        continue;
      }
      if(maxrate > 0.0 && rate > maxrate) {
        keep_tracklet[t] = false;
        n_rate_high++;
        continue;
      }
    }

    // --- minnights filter ---
    if(minnights > 0 && total_distinct_nights < minnights) {
      keep_tracklet[t] = false;
      n_nights++;
      continue;
    }

    // --- Per-detection filters (require looking up detection data via trk2det) ---
    // Collect detection indices for this tracklet
    long pos = t2d_offsets[t];
    // We need at least the first two detections for endpoint-based filters
    vector<long> det_indices;
    while(pos < (long)trk2det.size() && trk2det[pos].i1 == t) {
      det_indices.push_back(trk2det[pos].i2);
      pos++;
    }
    if(det_indices.size() < 2) {
      // Degenerate tracklet — skip silently
      keep_tracklet[t] = false;
      continue;
    }

    // For endpoint-based filters, use the first and last detection indices
    long d1 = det_indices[0];
    long d2 = det_indices[det_indices.size()-1];

    // Validate indices
    if(d1 < 0 || d1 >= detnum || d2 < 0 || d2 >= detnum) {
      cerr << "WARNING: tracklet " << t << " has out-of-range detection index, skipping\n";
      keep_tracklet[t] = false;
      continue;
    }

    // --- det_qual filter: reject tracklet if EITHER endpoint has det_qual < threshold ---
    if(mindetqual > 0) {
      if(detvec[d1].det_qual < mindetqual || detvec[d2].det_qual < mindetqual) {
        keep_tracklet[t] = false;
        n_detqual++;
        continue;
      }
    }

    // --- maxdmag filter ---
    if(maxdmag > 0.0) {
      double dmag = fabs(detvec[d1].mag - detvec[d2].mag);
      if(dmag > maxdmag) {
        keep_tracklet[t] = false;
        n_dmag++;
        continue;
      }
    }

    // --- maximdensity filter: reject if EITHER endpoint's image has too many detections ---
    if(maximdensity > 0) {
      long im1 = detvec[d1].image;
      long im2 = detvec[d2].image;
      if((im1 >= 0 && im1 < imnum && image_detcount[im1] > maximdensity) ||
         (im2 >= 0 && im2 < imnum && image_detcount[im2] > maximdensity)) {
        keep_tracklet[t] = false;
        n_imdensity++;
        continue;
      }
    }

    // --- excludeknown filter: reject if EITHER endpoint has known_obj > 0 ---
    if(excludeknown) {
      if(detvec[d1].known_obj > 0 || detvec[d2].known_obj > 0) {
        keep_tracklet[t] = false;
        n_known++;
        continue;
      }
    }

    // --- maxtraillen filter: reject if BOTH endpoints have trail_len > threshold ---
    if(maxtraillen > 0.0) {
      if(detvec[d1].trail_len > maxtraillen && detvec[d2].trail_len > maxtraillen) {
        keep_tracklet[t] = false;
        n_trail++;
        continue;
      }
    }
  }

  // Report filter statistics
  long n_kept = 0;
  for(long t = 0; t < pairnum; t++) if(keep_tracklet[t]) n_kept++;
  long n_rejected = pairnum - n_kept;
  cout << "\nFilter statistics:\n";
  cout << "  Input tracklets: " << pairnum << "\n";
  cout << "  Kept: " << n_kept << "\n";
  cout << "  Rejected: " << n_rejected << "\n";
  if(minrate > 0.0)   cout << "    Too slow (< " << minrate << " deg/day): " << n_rate_low << "\n";
  if(maxrate > 0.0)   cout << "    Too fast (> " << maxrate << " deg/day): " << n_rate_high << "\n";
  if(minnights > 0)   cout << "    Insufficient nights: " << n_nights << "\n";
  if(mindetqual > 0)  cout << "    Low det_qual (< " << mindetqual << "): " << n_detqual << "\n";
  if(maxdmag > 0.0)   cout << "    Large dmag (> " << maxdmag << " mag): " << n_dmag << "\n";
  if(maximdensity > 0) cout << "    Dense image (> " << maximdensity << " dets): " << n_imdensity << "\n";
  if(excludeknown)    cout << "    Known object: " << n_known << "\n";
  if(maxtraillen > 0.0) cout << "    Both endpoints trailed (> " << maxtraillen << " arcsec): " << n_trail << "\n";

  // ------------------------------------------------------------------
  // Pass 3: Build filtered output
  // ------------------------------------------------------------------

  // Collect all detection indices referenced by kept tracklets
  vector<bool> keep_detection(detnum, false);
  for(long t = 0; t < pairnum; t++) {
    if(!keep_tracklet[t]) continue;
    long pos = t2d_offsets[t];
    while(pos < (long)trk2det.size() && trk2det[pos].i1 == t) {
      long d = trk2det[pos].i2;
      if(d >= 0 && d < detnum) keep_detection[d] = true;
      pos++;
    }
  }

  // Build remapping from old detection index to new index
  vector<long> det_remap(detnum, -1);
  long new_detnum = 0;
  for(long d = 0; d < detnum; d++) {
    if(keep_detection[d]) {
      det_remap[d] = new_detnum++;
    }
  }

  // Build remapping from old tracklet index to new index
  vector<long> trk_remap(pairnum, -1);
  long new_trknum = 0;
  for(long t = 0; t < pairnum; t++) {
    if(keep_tracklet[t]) {
      trk_remap[t] = new_trknum++;
    }
  }

  cout << "\nOutput summary:\n";
  cout << "  Input detections: " << detnum << ", kept: " << new_detnum << "\n";
  cout << "  Input tracklets: " << pairnum << ", kept: " << new_trknum << "\n";

  // Write filtered pairdets
  ofstream outstream1;
  outstream1.open(outpairdetfile);
  if(!outstream1) {
    cerr << "ERROR: cannot open output pairdets file " << outpairdetfile << "\n";
    return(1);
  }
  outstream1 << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
  for(long d = 0; d < detnum; d++) {
    if(!keep_detection[d]) continue;
    outstream1 << fixed << setprecision(7) << detvec[d].MJD << "," << detvec[d].RA << "," << detvec[d].Dec << ",";
    outstream1 << fixed << setprecision(4) << detvec[d].mag << ",";
    outstream1 << fixed << setprecision(2) << detvec[d].trail_len << "," << detvec[d].trail_PA << ",";
    outstream1 << fixed << setprecision(4) << detvec[d].sigmag << ",";
    outstream1 << fixed << setprecision(3) << detvec[d].sig_across << "," << detvec[d].sig_along << ",";
    outstream1 << detvec[d].image << "," << detvec[d].idstring << "," << detvec[d].band << ",";
    outstream1 << detvec[d].obscode << "," << detvec[d].known_obj << ",";
    outstream1 << detvec[d].det_qual << "," << detvec[d].index << "\n";
  }
  outstream1.close();
  cout << "Wrote " << new_detnum << " detections to " << outpairdetfile << "\n";

  // Write filtered tracklets
  outstream1.open(outtrkfile);
  if(!outstream1) {
    cerr << "ERROR: cannot open output tracklets file " << outtrkfile << "\n";
    return(1);
  }
  outstream1 << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
  long out_trkct = 0;
  for(long t = 0; t < pairnum; t++) {
    if(!keep_tracklet[t]) continue;
    outstream1 << fixed << setprecision(7) << tracklets[t].Img1 << "," << tracklets[t].RA1 << "," << tracklets[t].Dec1 << ",";
    outstream1 << fixed << setprecision(7) << tracklets[t].Img2 << "," << tracklets[t].RA2 << "," << tracklets[t].Dec2 << ",";
    outstream1 << tracklets[t].npts << "," << out_trkct << "\n";
    out_trkct++;
  }
  outstream1.close();
  cout << "Wrote " << new_trknum << " tracklets to " << outtrkfile << "\n";

  // Write filtered trk2det (with remapped indices)
  outstream1.open(outtrk2detfile);
  if(!outstream1) {
    cerr << "ERROR: cannot open output trk2det file " << outtrk2detfile << "\n";
    return(1);
  }
  outstream1 << "#trk_ID,detnum\n";
  long out_t2d_count = 0;
  for(long t = 0; t < pairnum; t++) {
    if(!keep_tracklet[t]) continue;
    long new_trkid = trk_remap[t];
    long pos = t2d_offsets[t];
    while(pos < (long)trk2det.size() && trk2det[pos].i1 == t) {
      long d = trk2det[pos].i2;
      if(d >= 0 && d < detnum && det_remap[d] >= 0) {
        outstream1 << new_trkid << "," << det_remap[d] << "\n";
        out_t2d_count++;
      }
      pos++;
    }
  }
  outstream1.close();
  cout << "Wrote " << out_t2d_count << " trk2det entries to " << outtrk2detfile << "\n";

  return(0);
}
