// make_tracklets_omp.cpp
//
// OpenMP-parallel wrapper around make_tracklets7().
//
// Because find_pairs() only pairs images within [mintime, maxtime] (both
// sub-nightly for normal ATLAS settings), detections on different nights
// can never form a tracklet.  Nightly chunks of the detection catalog are
// therefore fully independent and can be processed in parallel.
//
// New argument vs. make_tracklets:
//   -nw N          number of OpenMP worker threads (default 1 = serial,
//                  behaviour identical to make_tracklets)
//   -noon_local H  snap night boundaries to local noon at UTC offset H hours
//                  (e.g. -10 for Hawaii HST).  Default: MJD integer floor.
//
// All other arguments are identical to make_tracklets.  The four output
// files (outimgs, pairdets, tracklets, trk2det) are written in exactly
// the same format.

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <algorithm>
#define NUMPOS 3

#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 5
#define BANDCOL 6
#define OBSCODECOL 7
#define COLS_TO_READ 14
#define MAXVEL 1.5
#define MAXTIME (1.5/24.0)
#define IMAGERAD 2.0
#define MAX_GCR 0.5
#define DEBUG 0
#define DEBUGA 0
#define DEBUGB 0

#define TESTOBSLON 203.74409
#define TESTPLXCOS 0.936241
#define TESTPLXSIN 0.351543
#define TESTOBSCODE "F51"
#define OBSCHECKTOL 1.0e-6

#define MINEARTHDIST 1.44e8
#define MAXEARTHDIST 1.55e8
#define MAXEARTHZ 1.0e5

static void show_usage()
{
  cerr << "Usage: make_tracklets_omp -dets detfile -imgs imfile -outimgs output image file/ \n";
  cerr << "-pairdets paired detection file -tracklets tracklet file -trk2det output tracklet-to-detection file/ \n";
  cerr << "-colformat column format file -imrad image radius(deg)/ \n";
  cerr << "-matchrad radius for image overlap calculation (deg)/ \n";
  cerr << "-trkfrac min fraction of overlapping images that must be included in a valid tracklet/ \n";
  cerr << "-maxtime max inter-image time interval (hr) -mintime min inter-image time interval (hr)/ \n";
  cerr << "-maxGCR maximum GRC -mintrkpts min. num. of tracklet points/ \n";
  cerr << "-max_netl maximum number of points for a non-exclusive (overlap permitted) tracklet/ \n";
  cerr << "-time_offset offset in seconds to be added to observations times to get UTC/ \n";
  cerr << "-minvel minimum angular velocity (deg/day) -maxvel maximum angular velocity (deg/day)/ \n";
  cerr << "-minarc minimum total angular arc (arcsec) -earth earthfile -obscode obscodefile -forcerun\n";
  cerr << "-nw num_omp_threads (default 1)\n";
  cerr << "-noon_local utc_offset_hours (snap night boundaries to local noon; default: MJD integer floor)\n";
  cerr << "\nor, at minimum\n\n";
  cerr << "make_tracklets_omp -dets detfile -earth earthfile -obscode obscodefile\n";
}


int main(int argc, char *argv[])
{
  vector <hldet> detvec = {};
  vector <observatory> observatory_list = {};
  vector <hlimage> img_log = {};

  vector <point3d> Earthpos;
  vector <point3d> Earthvel;
  vector <double> EarthMJD;
  string lnfromfile;
  int status = 0;
  long i = 0;
  int imct = 0;
  string indetfile;
  string inimfile;
  string earthfile;
  string obscodefile;
  string colformatfile;
  string outimfile = "outimfile01.txt";
  string pairdetfile = "pairdetfile01.csv";
  string trackletfile = "trackletfile01.csv";
  string trk2detfile = "trk2detfile01.csv";
  int idcol = IDCOL;
  int mjdcol = MJDCOL;
  int racol = RACOL;
  int deccol = DECCOL;
  int magcol = MAGCOL;
  int bandcol = BANDCOL;
  int obscodecol = OBSCODECOL;
  int trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col;
  trail_len_col = trail_PA_col = sigmag_col = sig_across_col = -1;
  sig_along_col = known_obj_col = det_qual_col = -1;
  int colreadct = 0;
  ifstream instream1;
  ofstream outstream1;
  string stest;
  int inimfile_set, colformatfile_set;
  inimfile_set = colformatfile_set = 0;
  int outimfile_default, pairdetfile_default, trackletfile_default, trk2detfile_default, imagerad_default;
  int maxtime_default, mintime_default, minvel_default, maxvel_default, matchrad_default, trkfrac_default;
  int maxgcr_default, minarc_default, mintrkpts_default, time_offset_default, maxnetl_default;
  MakeTrackletsConfig config;
  int n_workers = 1;
  double noon_local_offset = NAN;  // NAN means disabled; use MJD integer floor

  outimfile_default = pairdetfile_default = trackletfile_default = trk2detfile_default = imagerad_default = 1;
  maxtime_default = mintime_default = minvel_default = maxvel_default = matchrad_default = trkfrac_default = 1;
  maxgcr_default = minarc_default = mintrkpts_default = time_offset_default = maxnetl_default = 1;

  if(argc < 7) {
    show_usage();
    return(1);
  }

  i = 1;
  while(i < argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det" || string(argv[i]) == "--detection" || string(argv[i]) == "--detections") {
      if(i+1 < argc) {
        indetfile = argv[++i];
        i++;
      } else {
        cerr << "Detection file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img" || string(argv[i]) == "--image" || string(argv[i]) == "--images") {
      if(i+1 < argc) {
        inimfile = argv[++i];
        inimfile_set = 1;
        i++;
      } else {
        cerr << "Image file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-outim" || string(argv[i]) == "-outimg" || string(argv[i]) == "-outimgs" || string(argv[i]) == "--outimages" || string(argv[i]) == "--outimage" || string(argv[i]) == "--outimgs" || string(argv[i]) == "--outimg") {
      if(i+1 < argc) {
        outimfile = argv[++i];
        outimfile_default = 0;
        i++;
      } else {
        cerr << "Output image file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-pairdet" || string(argv[i]) == "-pairdets" || string(argv[i]) == "-detpairs") {
      if(i+1 < argc) {
        pairdetfile = argv[++i];
        pairdetfile_default = 0;
        i++;
      } else {
        cerr << "Output paired detection file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-tracklets" || string(argv[i]) == "-tracklet" || string(argv[i]) == "-pair" || string(argv[i]) == "-pairfile" || string(argv[i]) == "-pairs") {
      if(i+1 < argc) {
        trackletfile = argv[++i];
        trackletfile_default = 0;
        i++;
      } else {
        cerr << "Output pair file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-trk2det" || string(argv[i]) == "-t2d" || string(argv[i]) == "-tracklet2detection" || string(argv[i]) == "-trk2detfile" || string(argv[i]) == "--tracklet2detection") {
      if(i+1 < argc) {
        trk2detfile = argv[++i];
        trk2detfile_default = 0;
        i++;
      } else {
        cerr << "Output pair file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-imrad" || string(argv[i]) == "-imagerad") {
      if(i+1 < argc) {
        config.imagerad = stod(argv[++i]);
        imagerad_default = 0;
        i++;
        if(!isnormal(config.imagerad) || config.imagerad <= 0.0) {
          cerr << "Error: invalid image radius (" << config.imagerad << " deg) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Output image radius keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-matchrad" || string(argv[i]) == "--matchrad") {
      if(i+1 < argc) {
        config.matchrad = stod(argv[++i]);
        matchrad_default = 0;
        i++;
        if(!isnormal(config.matchrad) || config.matchrad <= 0.0) {
          cerr << "Error: invalid matching radius (" << config.matchrad << " deg) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Matching radius keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-trkfrac" || string(argv[i]) == "--trkfrac") {
      if(i+1 < argc) {
        config.trkfrac = stod(argv[++i]);
        trkfrac_default = 0;
        i++;
        if(!isnormal(config.trkfrac) || config.trkfrac <= 0.0) {
          cerr << "Error: invalid tracklet detection fraction (" << config.trkfrac << ") supplied.\n";
          return(2);
        }
      } else {
        cerr << "Tracklet detection fraction keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-maxtime") {
      if(i+1 < argc) {
        config.maxtime = stod(argv[++i]);
        maxtime_default = 0;
        i++;
        if(isnormal(config.maxtime) && config.maxtime > 0.0) {
          config.maxtime /= 24.0;
        } else {
          cerr << "Error: invalid maximum inter-image time interval (" << config.maxtime << " hr) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Maximum inter-image time interval keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-mintime") {
      if(i+1 < argc) {
        config.mintime = stod(argv[++i]);
        mintime_default = 0;
        i++;
        if((isnormal(config.mintime) || config.mintime == 0.0) && config.mintime >= 0.0) {
          config.mintime /= 24.0;
        } else {
          cerr << "Error: invalid minimum inter-image time interval (" << config.mintime << " hr) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Minimum inter-image time interval keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-minvel") {
      if(i+1 < argc) {
        config.minvel = stod(argv[++i]);
        minvel_default = 0;
        i++;
        if(!isnormal(config.minvel) && config.minvel != 0.0l) {
          cerr << "Error: invalid minimum angular velocity (" << config.minvel << " deg/day) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Minimum angular velocity keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-maxvel") {
      if(i+1 < argc) {
        config.maxvel = stod(argv[++i]);
        maxvel_default = 0;
        i++;
        if(!isnormal(config.maxvel) || config.maxvel <= 0.0) {
          cerr << "Error: invalid maximum angular velocity (" << config.maxvel << " deg/day) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Maximum angular velocity keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-maxGCR" || string(argv[i]) == "-maxgcr") {
      if(i+1 < argc) {
        config.maxgcr = stod(argv[++i]);
        maxgcr_default = 0;
        i++;
        if(!isnormal(config.maxgcr) || config.maxgcr <= 0.0) {
          cerr << "Error: invalid maximum Great Circle residual (" << config.maxgcr << " arcsec) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Output maximum Great Circle Residual keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-minarc") {
      if(i+1 < argc) {
        config.minarc = stod(argv[++i]);
        minarc_default = 0;
        i++;
        if(!isnormal(config.minarc) && config.minarc != 0.0l) {
          cerr << "Error: invalid minimum angular arc (" << config.minarc << " arcsec) supplied.\n";
          return(2);
        }
      } else {
        cerr << "Minimum angular arc keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earthfile" || string(argv[i]) == "--Earthfile" || string(argv[i]) == "--earth" || string(argv[i]) == "--Earth") {
      if(i+1 < argc) {
        earthfile = argv[++i];
        i++;
      } else {
        cerr << "Earth file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-mintrkpts" || string(argv[i]) == "-mpt" || string(argv[i]) == "-mintrackpts" || string(argv[i]) == "-minpts" || string(argv[i]) == "--minimumtrackletpoints" || string(argv[i]) == "--mintrackpoints" || string(argv[i]) == "--mintrackletpoints") {
      if(i+1 < argc) {
        config.mintrkpts = stoi(argv[++i]);
        mintrkpts_default = 0;
        i++;
      } else {
        cerr << "Min. tracklet points keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-max_netl" || string(argv[i]) == "-maxnetl" || string(argv[i]) == "-maxNETL" || string(argv[i]) == "-max_NETL" || string(argv[i]) == "--maximum_non-exclusive_tracklet_length" || string(argv[i]) == "--maximum_non_exclusive_tracklet_length" || string(argv[i]) == "--max_netl") {
      if(i+1 < argc) {
        config.max_netl = stoi(argv[++i]);
        maxnetl_default = 0;
        i++;
      } else {
        cerr << "Max. non-exclusive tracklet length keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-time_offset" || string(argv[i]) == "-offset" || string(argv[i]) == "-timeoffset" || string(argv[i]) == "-timeoff" || string(argv[i]) == "--time_offset" || string(argv[i]) == "--timeoffset" || string(argv[i]) == "--timeoff") {
      if(i+1 < argc) {
        config.time_offset = stod(argv[++i]);
        time_offset_default = 0;
        i++;
      } else {
        cerr << "Time offset keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" || string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" || string(argv[i]) == "--obscodes" || string(argv[i]) == "--observatorycodes") {
      if(i+1 < argc) {
        obscodefile = argv[++i];
        i++;
      } else {
        cerr << "Observatory code file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format" || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat" || string(argv[i]) == "--cformat") {
      if(i+1 < argc) {
        colformatfile = argv[++i];
        colformatfile_set = 1;
        i++;
      } else {
        cerr << "Column format file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-use_lowmem" || string(argv[i]) == "-lowmem" || string(argv[i]) == "-ulm" || string(argv[i]) == "-mem" || string(argv[i]) == "--use_lowmem" || string(argv[i]) == "--lowmem" || string(argv[i]) == "-incremental") {
      if(i+1 < argc) {
        config.use_lowmem = stoi(argv[++i]);
        i++;
      } else {
        cerr << "Keyword for choosing low-memory algorithm supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "-VERBOSE" || string(argv[i]) == "-VERB" || string(argv[i]) == "--verbose" || string(argv[i]) == "--VERBOSE" || string(argv[i]) == "--VERB") {
      if(i+1 < argc) {
        config.verbose = stoi(argv[++i]);
        i++;
      } else {
        cerr << "Verbosity keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-forcerun" || string(argv[i]) == "-force" || string(argv[i]) == "-fr" || string(argv[i]) == "-f" || string(argv[i]) == "--force" || string(argv[i]) == "--forcerun") {
      config.forcerun = 1;
      i++;
    } else if(string(argv[i]) == "-nw" || string(argv[i]) == "-n_workers" || string(argv[i]) == "-nworkers" || string(argv[i]) == "-nthreads" || string(argv[i]) == "--n_workers" || string(argv[i]) == "--nworkers" || string(argv[i]) == "--nthreads") {
      if(i+1 < argc) {
        n_workers = stoi(argv[++i]);
        i++;
      } else {
        cerr << "-nw keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-noon_local" || string(argv[i]) == "--noon_local") {
      if(i+1 < argc) {
        noon_local_offset = stod(argv[++i]);
        i++;
      } else {
        cerr << "-noon_local keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] << "\n";
      i++;
    }
  }

  if(indetfile.size() <= 0) {
    cerr << "Please supply an input detection file:\n\n";
    show_usage();
    return(1);
  }
  if(earthfile.size() <= 0) {
    cerr << "Please supply a heliocentric ephemeris file for the Earth:\n\n";
    show_usage();
    return(1);
  }
  if(obscodefile.size() <= 0) {
    cerr << "Please supply an observatory code file:\n\n";
    show_usage();
    return(1);
  }

  if(config.mintrkpts < 2) config.mintrkpts = 2;

  // Echo configuration
  cout << "\nInput detection file is called " << indetfile << "\n";
  if(inimfile_set == 1) cout << "Input image file = " << inimfile << "\n";
  else cout << "No input image file specified: image catalog will be generated internally.\n";
  if(outimfile_default == 0) cout << "Output image file will be called " << outimfile << "\n";
  else cout << "Defaulting to output image file name = " << outimfile << "\n";
  if(pairdetfile_default == 0) cout << "Output paired detection file will be called " << pairdetfile << "\n";
  else cout << "Defaulting to output paired detection file name = " << pairdetfile << "\n";
  if(trackletfile_default == 0) cout << "Output tracklet file will be called " << trackletfile << "\n";
  else cout << "Defaulting to output tracklet file name = " << trackletfile << "\n";
  if(trk2detfile_default == 0) cout << "Output tracklet-to-detection file will be called " << trk2detfile << "\n";
  else cout << "Defaulting to output tracklet-to-detection file name = " << trk2detfile << "\n";
  cout << "Heliocentric ephemeris file for the Earth: " << earthfile << "\n";

  if(n_workers > 1) {
    omp_set_num_threads(n_workers);
    cout << "OpenMP worker threads: " << n_workers << " (via -nw)\n";
  } else {
    n_workers = 1;
    cout << "Serial mode: n_workers=1 (use -nw N to parallelize)\n";
  }
  if(!isnan(noon_local_offset)) {
    cout << "Night boundaries: local noon at UTC offset " << noon_local_offset << " hours\n";
  } else {
    cout << "Night boundaries: MJD integer floor (default)\n";
  }

  // Read the column formatting file, if any
  if(colformatfile.size() > 0) {
    idcol = mjdcol = racol = deccol = magcol = bandcol = obscodecol = -1;
    instream1.open(colformatfile);
    if(!instream1) {
      cerr << "ERROR: unable to open input file " << colformatfile << "\n";
      return(1);
    }
    colreadct = 0;
    while(!instream1.eof() && !instream1.fail() && !instream1.bad() && colreadct < COLS_TO_READ) {
      instream1 >> stest;
      if(stest == "MJDCOL") { instream1 >> mjdcol; if(!instream1.fail()) colreadct++; }
      else if(stest == "RACOL") { instream1 >> racol; if(!instream1.fail()) colreadct++; }
      else if(stest == "DECCOL") { instream1 >> deccol; if(!instream1.fail()) colreadct++; }
      else if(stest == "MAGCOL") { instream1 >> magcol; if(!instream1.fail()) colreadct++; }
      else if(stest == "TRAILLENCOL") { instream1 >> trail_len_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "TRAILPACOL") { instream1 >> trail_PA_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "SIGMAGCOL") { instream1 >> sigmag_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "SIGACROSSCOL") { instream1 >> sig_across_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "SIGALONGCOL") { instream1 >> sig_along_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "IDCOL") { instream1 >> idcol; if(!instream1.fail()) colreadct++; }
      else if(stest == "BANDCOL") { instream1 >> bandcol; if(!instream1.fail()) colreadct++; }
      else if(stest == "OBSCODECOL") { instream1 >> obscodecol; if(!instream1.fail()) colreadct++; }
      else if(stest == "KNOWNOBJCOL") { instream1 >> known_obj_col; if(!instream1.fail()) colreadct++; }
      else if(stest == "DETQUALCOL") { instream1 >> det_qual_col; if(!instream1.fail()) colreadct++; }
      else { cout << "WARNING: unrecognized string " << stest << " read from column formatting file\n"; }
    }
    instream1.close();
    if(colreadct < COLS_TO_READ) {
      cout << "WARNING: only " << colreadct << " column specifications, of " << COLS_TO_READ << " expected, were read from column format file " << colformatfile << ".\n";
    }
  }

  // Read observatory code file
  status = read_obscode_file2(obscodefile, observatory_list, config.verbose);
  if(status != 0) {
    cerr << "ERROR reading observatory code file " << obscodefile << "\n";
    return(1);
  }
  if(config.forcerun <= 0) {
    char obscode[MINSTRINGLEN];
    double obslon, plxcos, plxsin;
    obslon = plxcos = plxsin;
    stringncopy01(obscode, TESTOBSCODE, MINSTRINGLEN);
    obscode_lookup(observatory_list, obscode, obslon, plxcos, plxsin);
    if(fabs(obslon-TESTOBSLON) > OBSCHECKTOL || fabs(plxcos-TESTPLXCOS) > OBSCHECKTOL || fabs(plxsin-TESTPLXSIN) > OBSCHECKTOL) {
      cerr << "ERROR: correct data from ObsCode F51 (Pan-STARRS 1) not read from file " << obscodefile << "\n";
      cerr << "If you are sure you want to run with the existing file, use -forcerun 1.\n";
      return(1);
    }
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file " << obscodefile << "\n";

  // Read input detection file
  status = read_detection_filemt2(indetfile, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol, trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col, detvec, config.verbose, config.forcerun);
  if(status == 0) {
    cout << "Input file " << indetfile << " read successfully.\n";
  } else if(status == 1) {
    cerr << "Warning: file read failed\n";
    return(1);
  } else if(status == 2) {
    cerr << "Warning: file possibly corrupted\n";
    return(2);
  } else {
    cerr << "Warning: unknown file read problem\n";
    return(3);
  }

  // Time-sort the detection vector
  sort(detvec.begin(), detvec.end(), early_hldet());

  // Get image information, if there is an image file
  if(inimfile.size() > 0) {
    cout << "About to read image file " << inimfile << "\n";
    status = read_image_file(inimfile, img_log);
    cout << "Read image file with " << img_log.size() << " lines\n";
    if(status == 0) {
      sort(img_log.begin(), img_log.end(), early_hlimage());
    } else {
      cerr << "Warning: failed to read supplied image file " << inimfile << "\n";
      cerr << "Constructing image table by inference from input detections instead\n";
      img_log = {};
    }
  }

  EarthMJD = {};
  Earthpos = {};
  Earthvel = {};
  read_horizons_csv(earthfile, EarthMJD, Earthpos, Earthvel);
  cout << "Finished reading heliocentric ephemeris file " << earthfile << " for the Earth.\n";
  for(i = 0; i < long(Earthpos.size()); i++) {
    double rho = sqrt(Earthpos[i].x*Earthpos[i].x + Earthpos[i].y*Earthpos[i].y);
    if(rho < MINEARTHDIST || rho > MAXEARTHDIST || fabs(Earthpos[i].z) > MAXEARTHZ) {
      cerr << "ERROR: Impossible earth coordinates from file " << earthfile << "\n";
      return(2);
    }
  }

  status = load_image_table(img_log, detvec, config.time_offset, observatory_list, EarthMJD, Earthpos, Earthvel);

  // Populate detvec[i].image from the freshly-built img_log so that the
  // nightly partitioning loop below can assign detections to nights.
  // We use img_log[j].startind/endind (set by load_image_table) directly
  // rather than calling load_image_indices(), because load_image_indices()
  // uses a strict MJD-tolerance re-match that can fail on multi-observatory
  // edge cases where load_image_table computed a mean image MJD that drifts
  // slightly outside the per-detection tolerance.  The startind/endind-based
  // assignment is always consistent with the image grouping actually used.
  for(long ji = 0; ji < long(img_log.size()); ji++) {
    for(long di = img_log[ji].startind; di < img_log[ji].endind; di++) {
      if(di >= 0 && di < long(detvec.size())) {
        detvec[di].image = ji;
      }
    }
  }
  {
    long n_indexed = 0;
    for(long di = 0; di < long(detvec.size()); di++) if(detvec[di].image >= 0) n_indexed++;
    cout << "Assigned image indices: " << n_indexed << " of " << detvec.size()
         << " detections indexed across " << img_log.size() << " images.\n";
  }

  // Replace any invalid exposure times with the default value
  long exp_resetnum = 0;
  for(i = 0; i < long(img_log.size()); i++) {
    if(img_log[i].exptime <= 0.0l) {
      if(config.verbose > 0) cout << "Correcting exposure time on image " << i << "\n";
      img_log[i].exptime = config.exptime;
      exp_resetnum++;
    }
  }
  cout << "In main, exposure time was corrected for " << exp_resetnum << " out of " << img_log.size() << " images\n";

  // Write image log (same as make_tracklets, done before the parallel section)
  cout << "Writing output image catalog " << outimfile << " with " << img_log.size() << " lines\n";
  outstream1.open(outimfile);
  for(imct = 0; imct < long(img_log.size()); imct++) {
    outstream1 << fixed << setprecision(8) << img_log[imct].MJD << " " << img_log[imct].RA;
    outstream1 << fixed << setprecision(8) << " " << img_log[imct].Dec << " " << img_log[imct].obscode << " ";
    outstream1 << fixed << setprecision(1) << img_log[imct].X << " " << img_log[imct].Y << " " << img_log[imct].Z << " ";
    outstream1 << fixed << setprecision(4) << img_log[imct].VX << " " << img_log[imct].VY << " " << img_log[imct].VZ << " ";
    outstream1 << img_log[imct].startind << " " << img_log[imct].endind << " " << img_log[imct].exptime << "\n";
  }
  outstream1.close();

  // -----------------------------------------------------------------------
  // Nightly partition and parallel make_tracklets7 calls
  // -----------------------------------------------------------------------

  long nimages = img_log.size();
  long ndet = detvec.size();

  if(nimages == 0 || ndet == 0) {
    cout << "WARNING: no images or no detections; writing empty output files.\n";
    // Fall through to write empty output files.
  }

  // Compute night boundaries.
  // Strategy: assign each image an integer "night index" based on MJD.
  // If noon_local_offset is set, snap boundaries to local noon (matching
  // split_tracklets_by_time logic).  Otherwise use MJD integer floor.
  double win_origin = 0.0;
  if(!isnan(noon_local_offset) && nimages > 0) {
    double noon_utc_frac = fmod((12.0 - noon_local_offset) / 24.0, 1.0);
    if(noon_utc_frac < 0.0) noon_utc_frac += 1.0;
    double min_mjd = img_log[0].MJD;
    for(long ii = 1; ii < nimages; ii++) {
      if(img_log[ii].MJD < min_mjd) min_mjd = img_log[ii].MJD;
    }
    win_origin = floor(min_mjd) + noon_utc_frac;
    if(win_origin > min_mjd) win_origin -= 1.0;
    cout << "Noon-local snap: UTC offset=" << noon_local_offset
         << " h, window origin MJD=" << fixed << setprecision(5) << win_origin << "\n";
  } else {
    // Use MJD integer floor: win_origin = 0 means night_index = floor(MJD)
    win_origin = NAN;
  }

  // Assign each image a night index
  vector<long> img_night(nimages);
  for(long ii = 0; ii < nimages; ii++) {
    if(!isnan(win_origin)) {
      img_night[ii] = (long)floor(img_log[ii].MJD - win_origin);
    } else {
      img_night[ii] = (long)floor(img_log[ii].MJD);
    }
  }

  // Collect unique night indices (sorted)
  vector<long> night_ids = img_night;
  sort(night_ids.begin(), night_ids.end());
  night_ids.erase(unique(night_ids.begin(), night_ids.end()), night_ids.end());
  long nnights = night_ids.size();
  cout << "Found " << nnights << " distinct nights; running with " << n_workers << " threads.\n";

  // For each night, collect the image indices
  vector<vector<long>> night_imgs(nnights);
  for(long ii = 0; ii < nimages; ii++) {
    long nidx = (long)(lower_bound(night_ids.begin(), night_ids.end(), img_night[ii]) - night_ids.begin());
    night_imgs[nidx].push_back(ii);
  }

  // Per-night output containers (indexed by night)
  vector<vector<hldet>>    night_pairdets(nnights);
  vector<vector<tracklet>> night_tracklets(nnights);
  vector<vector<longpair>> night_trk2det(nnights);
  vector<int>              night_status(nnights, 0);

  // Parallel loop over nights
  #pragma omp parallel for num_threads(n_workers) schedule(dynamic)
  for(long ni = 0; ni < nnights; ni++) {
    // Build per-night image sub-log
    vector<hlimage> img_log_night;
    img_log_night.reserve(night_imgs[ni].size());
    for(long old_img : night_imgs[ni]) {
      img_log_night.push_back(img_log[old_img]);
    }

    // Build per-night detection sub-vector by iterating the per-night image
    // sub-log in order and collecting detections via each image's startind/endind
    // range in the global detvec.  This preserves the exact image-grouping order
    // established by load_image_table() and makes the per-night detvec order
    // consistent with img_log_night without a secondary sort.  A secondary sort
    // (early_hldet) is not safe here because early_hldet is not a strict weak
    // ordering when detections from different observatories have MJDs within
    // IMAGETIMETOL of each other; std::sort with such a comparator has undefined
    // behaviour and can produce orderings that are inconsistent with img_log.

    vector<hldet> detvec_night;
    long night_det_offset = 0;
    for(long ni_img = 0; ni_img < (long)night_imgs[ni].size(); ni_img++) {
      long old_img = night_imgs[ni][ni_img];
      long s = img_log[old_img].startind;
      long e = img_log[old_img].endind;
      // Update img_log_night[ni_img] startind/endind to per-night offsets.
      if(s < 0 || e <= s) {
        // Image has no detections (startind=endind=0 or invalid): keep as-is,
        // reset to 0 so load_image_indices inside make_tracklets7 handles it.
        img_log_night[ni_img].startind = night_det_offset;
        img_log_night[ni_img].endind   = night_det_offset;
      } else {
        img_log_night[ni_img].startind = night_det_offset;
        for(long di = s; di < e && di < ndet; di++) {
          hldet d = detvec[di];
          d.image = ni_img;
          detvec_night.push_back(d);
        }
        night_det_offset = (long)detvec_night.size();
        img_log_night[ni_img].endind = night_det_offset;
      }
    }

    if(detvec_night.empty()) {
      night_status[ni] = 0;
      continue;
    }

    // Call make_tracklets7 for this night
    vector<hldet>    pd_night;
    vector<tracklet> trk_night;
    vector<longpair> t2d_night;

    int st = make_tracklets7(detvec_night, img_log_night, config, pd_night, trk_night, t2d_night);
    night_status[ni] = st;
    if(st != 0) {
      cerr << "ERROR: make_tracklets7 failed with status " << st << " on night index " << ni << "\n";
      continue;
    }

    night_pairdets[ni]  = move(pd_night);
    night_tracklets[ni] = move(trk_night);
    night_trk2det[ni]   = move(t2d_night);
  }
  // End of parallel region

  // Check for any per-night failures
  for(long ni = 0; ni < nnights; ni++) {
    if(night_status[ni] != 0) {
      cerr << "ERROR: night " << ni << " failed; aborting merge.\n";
      return(night_status[ni]);
    }
  }

  // -----------------------------------------------------------------------
  // Serial offset pass: compute cumulative per-night pairdet / tracklet
  // start offsets.  O(nights), negligible cost.
  //
  // Index layout:
  //   pairdets:  global index = pd_offset[ni]  + night-local pd index
  //   tracklets: global index = trk_offset[ni] + night-local trk index
  //   trk2det entries: (global trk index, global pd index)
  // -----------------------------------------------------------------------
  vector<long> pd_offset(nnights, 0);
  vector<long> trk_offset(nnights, 0);
  {
    long pd_cumul = 0;
    long trk_cumul = 0;
    for(long ni = 0; ni < nnights; ni++) {
      pd_offset[ni]  = pd_cumul;
      trk_offset[ni] = trk_cumul;
      pd_cumul  += (long)night_pairdets[ni].size();
      trk_cumul += (long)night_tracklets[ni].size();
    }
    cout << "Offset pass complete: " << pd_cumul << " total pairdets, "
         << trk_cumul << " total tracklets from " << nnights << " nights.\n";
  }

  // -----------------------------------------------------------------------
  // Parallel apply-offsets + per-night write phase.
  //
  // Each worker:
  //   1. Applies global offsets to its in-memory vectors in-place.
  //   2. Writes per-night part files: <file>.part_<ni>
  //      (header row only on ni == 0; subsequent parts are headerless so
  //       the concat phase can simply stream them together without stripping).
  //   3. Clears its per-night vectors to release memory.
  // -----------------------------------------------------------------------
  cout << "Output image catalog " << outimfile << ", with " << img_log.size() << " lines, has been written\n";

  // We write into the same directory as the output files.  Part file naming:
  //   <pairdetfile>.part_<ni>
  //   <trackletfile>.part_<ni>
  //   <trk2detfile>.part_<ni>

  vector<int> write_status(nnights, 0);

  #pragma omp parallel for num_threads(n_workers) schedule(dynamic)
  for(long ni = 0; ni < nnights; ni++) {
    bool write_header = (ni == 0);

    // --- apply offsets to pairdets ---
    long poff = pd_offset[ni];
    long toff = trk_offset[ni];
    for(auto &pd : night_pairdets[ni]) {
      long local_img = pd.image;
      if(local_img >= 0 && local_img < (long)night_imgs[ni].size()) {
        pd.image = night_imgs[ni][local_img];
      }
      // pd.index is the original detection index, no remapping needed.
    }

    // --- apply offsets to tracklets ---
    for(long ti = 0; ti < (long)night_tracklets[ni].size(); ti++) {
      tracklet &trk = night_tracklets[ni][ti];
      long local_img1 = trk.Img1;
      long local_img2 = trk.Img2;
      if(local_img1 >= 0 && local_img1 < (long)night_imgs[ni].size()) {
        trk.Img1 = night_imgs[ni][local_img1];
      }
      if(local_img2 >= 0 && local_img2 < (long)night_imgs[ni].size()) {
        trk.Img2 = night_imgs[ni][local_img2];
      }
      trk.trk_ID = toff + ti;
    }

    // --- apply offsets to trk2det ---
    for(auto &tp : night_trk2det[ni]) {
      tp.i1 += toff;   // tracklet global index
      tp.i2 += poff;   // pairdet global index
    }

    // --- write per-night part files ---
    string pd_part   = pairdetfile   + ".part_" + to_string(ni);
    string trk_part  = trackletfile  + ".part_" + to_string(ni);
    string t2d_part  = trk2detfile   + ".part_" + to_string(ni);

    {
      ofstream os(pd_part);
      if(!os) { write_status[ni] = 1; goto cleanup_ni; }
      if(write_header)
        os << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,origindex\n";
      for(const auto &pd : night_pairdets[ni]) {
        os << fixed << setprecision(7) << pd.MJD << "," << pd.RA << "," << pd.Dec << ",";
        os << fixed << setprecision(4) << pd.mag << ",";
        os << fixed << setprecision(2) << pd.trail_len << "," << pd.trail_PA << ",";
        os << fixed << setprecision(4) << pd.sigmag << ",";
        os << fixed << setprecision(3) << pd.sig_across << "," << pd.sig_along << ",";
        os << pd.image << "," << pd.idstring << "," << pd.band << ",";
        os << pd.obscode << "," << pd.known_obj << ",";
        os << pd.det_qual << "," << pd.index << "\n";
      }
    }

    {
      ofstream os(trk_part);
      if(!os) { write_status[ni] = 1; goto cleanup_ni; }
      if(write_header)
        os << "#Image1,RA1,Dec1,Image2,RA2,Dec2,npts,trk_ID\n";
      for(const auto &trk : night_tracklets[ni]) {
        os << fixed << setprecision(7) << trk.Img1 << "," << trk.RA1 << "," << trk.Dec1 << ",";
        os << fixed << setprecision(7) << trk.Img2 << "," << trk.RA2 << "," << trk.Dec2 << ",";
        os << trk.npts << "," << trk.trk_ID << "\n";
      }
    }

    {
      ofstream os(t2d_part);
      if(!os) { write_status[ni] = 1; goto cleanup_ni; }
      if(write_header)
        os << "#trk_ID,detnum\n";
      for(const auto &tp : night_trk2det[ni]) {
        os << tp.i1 << "," << tp.i2 << "\n";
      }
    }

    cleanup_ni:
    // Release per-night memory immediately after writing
    { vector<hldet>().swap(night_pairdets[ni]); }
    { vector<tracklet>().swap(night_tracklets[ni]); }
    { vector<longpair>().swap(night_trk2det[ni]); }
  }
  // End of parallel write region

  for(long ni = 0; ni < nnights; ni++) {
    if(write_status[ni] != 0) {
      cerr << "ERROR: failed to write part files for night " << ni << "\n";
      return(1);
    }
  }

  // -----------------------------------------------------------------------
  // Serial concat phase: stream part files into final output files, then
  // delete the parts.  Done in C++ with a fixed-size read buffer to avoid
  // a shell dependency.
  // -----------------------------------------------------------------------
  auto concat_parts = [&](const string &outfile, const string &stem, long n) -> int {
    // open destination
    ofstream out(outfile, ios::binary);
    if(!out) {
      cerr << "ERROR: cannot open output file " << outfile << " for writing\n";
      return(1);
    }
    const size_t BUF = 1 << 22;  // 4 MiB read buffer
    vector<char> buf(BUF);
    for(long ni = 0; ni < n; ni++) {
      string part = stem + ".part_" + to_string(ni);
      ifstream in(part, ios::binary);
      if(!in) {
        cerr << "ERROR: cannot open part file " << part << " for reading\n";
        return(1);
      }
      while(in) {
        in.read(buf.data(), BUF);
        streamsize got = in.gcount();
        if(got > 0) out.write(buf.data(), got);
      }
      in.close();
      if(remove(part.c_str()) != 0) {
        cerr << "WARNING: could not remove part file " << part << "\n";
      }
    }
    out.close();
    return(0);
  };

  cout << "Concatenating per-night part files into final output files...\n";

  cout << "  Writing paired detection file " << pairdetfile << "\n";
  if(concat_parts(pairdetfile, pairdetfile, nnights) != 0) return(1);

  cout << "  Writing tracklet file " << trackletfile << "\n";
  if(concat_parts(trackletfile, trackletfile, nnights) != 0) return(1);

  cout << "  Writing trk2det file " << trk2detfile << "\n";
  if(concat_parts(trk2detfile, trk2detfile, nnights) != 0) return(1);

  cout << "All output files written successfully.\n";

  return(0);
}
