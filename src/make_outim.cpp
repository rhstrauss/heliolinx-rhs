// make_outim.cpp:
// Generates only the output image (outim) file from a detection catalog.
// This is a stripped-down version of make_tracklets that stops after writing
// the outim file, skipping all tracklet generation. Intended to prepare
// the outim file for Sorcha labelling without the expense of full tracklet creation.

#include "solarsyst_dyn_geo01.h"
#include "cmath"

#define IDCOL 1
#define MJDCOL 2
#define RACOL 3
#define DECCOL 4
#define MAGCOL 5
#define BANDCOL 6
#define OBSCODECOL 7
#define COLS_TO_READ 14

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
  cerr << "Usage: make_outim -dets detfile -outimgs output_image_file\n";
  cerr << "-earth earthfile -obscode obscodefile\n";
  cerr << "Optional: -imgs input_image_file -colformat column_format_file\n";
  cerr << "-time_offset offset_seconds -verbose verbosity -forcerun\n";
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
  int inimfile_set = 0;
  int colformatfile_set = 0;
  int outimfile_default = 1;
  MakeTrackletsConfig config;

  if(argc < 7) {
    show_usage();
    return(1);
  }

  i = 1;
  while(i < argc) {
    cout << "Checking out argv[" << i << "] = " << argv[i] << ".\n";
    if(string(argv[i]) == "-d" || string(argv[i]) == "-dets" || string(argv[i]) == "-det" || string(argv[i]) == "--dets" || string(argv[i]) == "--det") {
      if(i+1 < argc) {
        indetfile = argv[++i];
        i++;
      } else {
        cerr << "Detection file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-i" || string(argv[i]) == "-imgs" || string(argv[i]) == "-inimgs" || string(argv[i]) == "-img" || string(argv[i]) == "--inimgs" || string(argv[i]) == "--img") {
      if(i+1 < argc) {
        inimfile = argv[++i];
        inimfile_set = 1;
        i++;
      } else {
        cerr << "Image file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-outim" || string(argv[i]) == "-outimg" || string(argv[i]) == "-outimgs" || string(argv[i]) == "--outimages" || string(argv[i]) == "--outimg") {
      if(i+1 < argc) {
        outimfile = argv[++i];
        outimfile_default = 0;
        i++;
      } else {
        cerr << "Output image file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-earth" || string(argv[i]) == "-e" || string(argv[i]) == "-Earth" || string(argv[i]) == "--earth" || string(argv[i]) == "--earthfile") {
      if(i+1 < argc) {
        earthfile = argv[++i];
        i++;
      } else {
        cerr << "Earth file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-obscode" || string(argv[i]) == "-obs" || string(argv[i]) == "-oc" || string(argv[i]) == "-obscodes" || string(argv[i]) == "--obscode" || string(argv[i]) == "--obscodes") {
      if(i+1 < argc) {
        obscodefile = argv[++i];
        i++;
      } else {
        cerr << "Observatory code file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-colformat" || string(argv[i]) == "-format" || string(argv[i]) == "-col" || string(argv[i]) == "-cf" || string(argv[i]) == "-colfmt" || string(argv[i]) == "--colformat" || string(argv[i]) == "--columnformat") {
      if(i+1 < argc) {
        colformatfile = argv[++i];
        colformatfile_set = 1;
        i++;
      } else {
        cerr << "Column format file keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-time_offset" || string(argv[i]) == "-offset" || string(argv[i]) == "-timeoffset" || string(argv[i]) == "--time_offset" || string(argv[i]) == "--timeoffset") {
      if(i+1 < argc) {
        config.time_offset = stod(argv[++i]);
        i++;
      } else {
        cerr << "Time offset keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb" || string(argv[i]) == "--verbose") {
      if(i+1 < argc) {
        config.verbose = stoi(argv[++i]);
        i++;
      } else {
        cerr << "Verbosity keyword supplied with no corresponding argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-forcerun" || string(argv[i]) == "-force" || string(argv[i]) == "-f" || string(argv[i]) == "--forcerun" || string(argv[i]) == "--force") {
      config.forcerun = 1;
      i++;
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

  cout << "\nInput detection file: " << indetfile << "\n";
  if(inimfile_set == 1) cout << "Input image file: " << inimfile << "\n";
  else cout << "No input image file specified: image catalog will be generated from detections.\n";
  if(colformatfile_set == 1) cout << "Column formatting file: " << colformatfile << "\n";
  if(outimfile_default == 0) cout << "Output image file: " << outimfile << "\n";
  else cout << "Defaulting to output image file name: " << outimfile << "\n";
  cout << "Earth ephemeris file: " << earthfile << "\n";
  cout << "Observatory code file: " << obscodefile << "\n";

  // Read column format file
  if(colformatfile_set == 1) {
    idcol = mjdcol = racol = deccol = magcol = bandcol = obscodecol = -1;
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
      cout << "WARNING: only " << colreadct << " column specifications of " << COLS_TO_READ << " expected were read from " << colformatfile << ".\n";
    }
  }

  cout << "Column specifications: IDCOL=" << idcol << " MJDCOL=" << mjdcol
       << " RACOL=" << racol << " DECCOL=" << deccol
       << " MAGCOL=" << magcol << " BANDCOL=" << bandcol
       << " OBSCODECOL=" << obscodecol << "\n";

  // Read observatory codes
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
      cerr << "ERROR: correct data for ObsCode F51 (Pan-STARRS 1) not found in " << obscodefile << "\n";
      cerr << "Use -forcerun to override this check.\n";
      return(1);
    }
  }
  cout << "Read " << observatory_list.size() << " lines from observatory code file.\n";

  // Read detections
  status = read_detection_filemt2(indetfile, mjdcol, racol, deccol, magcol, idcol, bandcol, obscodecol, trail_len_col, trail_PA_col, sigmag_col, sig_across_col, sig_along_col, known_obj_col, det_qual_col, detvec, config.verbose, config.forcerun);
  if(status == 0) {
    cout << "Detection file " << indetfile << " read successfully.\n";
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

  // Time-sort detections
  sort(detvec.begin(), detvec.end(), early_hldet());

  // Read input image file if provided
  if(inimfile.size() > 0) {
    cout << "Reading input image file " << inimfile << "\n";
    status = read_image_file(inimfile, img_log);
    if(status == 0) {
      sort(img_log.begin(), img_log.end(), early_hlimage());
    } else {
      cerr << "Warning: failed to read image file " << inimfile << "; constructing from detections.\n";
      img_log = {};
    }
  }

  // Read Earth ephemeris
  EarthMJD = {};
  Earthpos = {};
  Earthvel = {};
  read_horizons_csv(earthfile, EarthMJD, Earthpos, Earthvel);
  cout << "Read Earth ephemeris from " << earthfile << ".\n";
  for(i = 0; i < long(Earthpos.size()); i++) {
    double rho = sqrt(Earthpos[i].x*Earthpos[i].x + Earthpos[i].y*Earthpos[i].y);
    if(rho < MINEARTHDIST || rho > MAXEARTHDIST || fabs(Earthpos[i].z) > MAXEARTHZ) {
      cerr << "ERROR: impossible Earth coordinates at index " << i << ": "
           << Earthpos[i].x << ", " << Earthpos[i].y << ", " << Earthpos[i].z << "\n";
      return(2);
    }
  }

  // Build image table with Earth positions
  status = load_image_table(img_log, detvec, config.time_offset, observatory_list, EarthMJD, Earthpos, Earthvel);

  // Fix any zero/negative exposure times
  long exp_resetnum = 0;
  for(i = 0; i < long(img_log.size()); i++) {
    if(img_log[i].exptime <= 0.0l) {
      img_log[i].exptime = config.exptime;
      exp_resetnum++;
    }
  }
  cout << "Exposure time corrected for " << exp_resetnum << " of " << img_log.size() << " images.\n";

  // Write outim file
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
  cout << "Output image catalog " << outimfile << " written successfully.\n";

  return(0);
}
