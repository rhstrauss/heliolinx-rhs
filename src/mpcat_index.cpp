// mpcat_index.cpp  (June 04, 2026)
//
// Build the pre-indexed, time-sorted binary catalog (mpcat.bin) consumed by
// mpcat_check. Streams one or more MPC 80-column observation files (e.g. the
// MPC bulk archives NumObs.txt and UnnObs.txt), parses every optical RA/Dec
// observation into a compact mpcdet record, sorts the whole set by MJD, and
// writes it as a raw mpcdet[] plus a small text sidecar (mpcat.bin.hdr) for
// validation and provenance.
//
// This is the heliolinx-native, build-once index that lets mpcat_check answer
// "which of these measurements already exist in the MPC?" by mmap'ing the
// catalog and binary-searching the relevant MJD window -- without ever loading
// the ~40 GB / ~539M-line archive into RAM at query time.
//
// Usage:
//   mpcat_index -numobs NumObs.txt -unnobs UnnObs.txt [-extra CmtObs.txt ...] \
//               -out mpcat.bin [-reserve 600000000]
//
// Any of -numobs / -unnobs / -extra may be omitted; at least one input file is
// required. -reserve pre-allocates the record vector to avoid reallocation
// spikes (default 600M ~ 34 GB for ~539M records). The build is memory-heavy
// (whole catalog held + sorted in RAM, ~30 GB) -- run on a big-RAM node.

#include "solarsyst_dyn_geo01.h"

#define WRITE_CHUNK 16000000L  // records per binary write block

static void show_usage()
{
  cerr << "Usage: mpcat_index -numobs NumObs.txt -unnobs UnnObs.txt \\\n";
  cerr << "       [-extra file1 [file2 ...]] -out mpcat.bin [-reserve N]\n";
  cerr << "\n";
  cerr << "Builds the time-sorted binary MPC catalog used by mpcat_check.\n";
  cerr << "\n";
  cerr << "Arguments:\n";
  cerr << "  -numobs   MPC 80-column file of numbered-object observations\n";
  cerr << "  -unnobs   MPC 80-column file of unnumbered-object observations\n";
  cerr << "  -extra    One or more additional MPC 80-column files (e.g. CmtObs.txt)\n";
  cerr << "  -out      Output binary catalog (a sidecar <out>.hdr is also written)\n";
  cerr << "  -reserve  Records to pre-allocate (default 600000000)\n";
  cerr << "  (At least one of -numobs / -unnobs / -extra must be supplied.)\n";
}

// Stream one MPC 80-column file, appending parsed optical observations to recs.
// Returns 0 on a clean read to EOF; reports counts of kept and skipped lines.
static int ingest_obs80(const string &infile, vector<mpcdet> &recs, long &kept, long &skipped)
{
  ifstream instream1(infile);
  if(!instream1) {
    cerr << "ERROR: can't open input file " << infile << "\n";
    return(1);
  }
  string lnfromfile;
  mpcdet rec;
  long lct = 0;
  kept = skipped = 0;
  while(getline(instream1, lnfromfile)) {
    lct++;
    if(mpcat_parse_obs80line(lnfromfile, rec) == 0) {
      recs.push_back(rec);
      kept++;
    } else {
      skipped++;
    }
    if(lct % 10000000L == 0) {
      cout << "  " << infile << ": " << lct/1000000L << "M lines ("
           << kept << " kept, " << skipped << " skipped)\n";
      cout.flush();
    }
  }
  cout << "Finished " << infile << ": " << lct << " lines, "
       << kept << " kept, " << skipped << " skipped\n";
  return(0);
}

int main(int argc, char *argv[])
{
  string numobs_file, unnobs_file, outfile;
  vector<string> extra_files = {};
  vector<string> source_files = {};
  long reserve_n = 600000000L;
  long i;

  if(argc < 3) { show_usage(); return(1); }

  i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-numobs" || string(argv[i]) == "--numobs") {
      if(i+1 < argc) { numobs_file = argv[++i]; i++; }
      else { cerr << "ERROR: -numobs requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-unnobs" || string(argv[i]) == "--unnobs") {
      if(i+1 < argc) { unnobs_file = argv[++i]; i++; }
      else { cerr << "ERROR: -unnobs requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-extra" || string(argv[i]) == "-in" || string(argv[i]) == "--extra") {
      i++;
      while(i < argc && argv[i][0] != '-') { extra_files.push_back(string(argv[i])); i++; }
      if(extra_files.empty()) { cerr << "ERROR: -extra requires at least one file argument\n"; return(1); }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "-o" || string(argv[i]) == "--out") {
      if(i+1 < argc) { outfile = argv[++i]; i++; }
      else { cerr << "ERROR: -out requires an argument\n"; return(1); }
    } else if(string(argv[i]) == "-reserve" || string(argv[i]) == "--reserve") {
      if(i+1 < argc) { reserve_n = stol(argv[++i]); i++; }
      else { cerr << "ERROR: -reserve requires an argument\n"; return(1); }
    } else {
      cerr << "Warning: unrecognized keyword " << argv[i] << "\n";
      i++;
    }
  }

  // Assemble the ordered list of source files.
  if(!numobs_file.empty()) source_files.push_back(numobs_file);
  if(!unnobs_file.empty()) source_files.push_back(unnobs_file);
  for(i = 0; i < long(extra_files.size()); i++) source_files.push_back(extra_files[i]);

  if(source_files.empty()) {
    cerr << "ERROR: no input files specified\n";
    show_usage();
    return(1);
  }
  if(outfile.empty()) {
    cerr << "ERROR: no output file specified (-out)\n";
    show_usage();
    return(1);
  }

  cout << "mpcat_index: building time-sorted MPC catalog\n";
  cout << "Record size: " << sizeof(mpcdet) << " bytes\n";
  cout << "Input files (" << source_files.size() << "):\n";
  for(i = 0; i < long(source_files.size()); i++) cout << "  " << source_files[i] << "\n";
  cout << "Output: " << outfile << " (+ " << outfile << ".hdr)\n";

  // Pre-allocate to avoid reallocation spikes during ingest.
  vector<mpcdet> recs;
  try {
    recs.reserve(reserve_n);
  } catch(const std::exception &ex) {
    cerr << "WARNING: could not reserve " << reserve_n << " records (" << ex.what()
         << "); continuing without reservation.\n";
  }

  long total_kept = 0, total_skipped = 0;
  for(i = 0; i < long(source_files.size()); i++) {
    long kept = 0, skipped = 0;
    int status = ingest_obs80(source_files[i], recs, kept, skipped);
    if(status != 0) {
      cerr << "ERROR ingesting " << source_files[i] << " (status " << status << ")\n";
      return(2);
    }
    total_kept += kept;
    total_skipped += skipped;
  }

  long nrec = long(recs.size());
  cout << "\nTotal records parsed: " << nrec << " (" << total_skipped << " lines skipped)\n";
  if(nrec == 0) {
    cerr << "ERROR: no records parsed; nothing to write.\n";
    return(2);
  }

  // Sort by MJD (ascending). This is the index: a contiguous, time-ordered
  // array that mpcat_check binary-searches for a query window.
  cout << "Sorting " << nrec << " records by MJD...\n";
  cout.flush();
  std::sort(recs.begin(), recs.end(),
            [](const mpcdet &a, const mpcdet &b) { return a.MJD < b.MJD; });
  double tmin = recs.front().MJD;
  double tmax = recs.back().MJD;
  cout << "Sorted. MJD range: " << fixed << setprecision(6) << tmin << " to " << tmax << "\n";

  // Write the binary catalog in chunks (robust for multi-GB output).
  ofstream outstream1(outfile, ios::binary);
  if(!outstream1) {
    cerr << "ERROR: can't open output file " << outfile << "\n";
    return(1);
  }
  long written = 0;
  while(written < nrec) {
    long block = nrec - written;
    if(block > WRITE_CHUNK) block = WRITE_CHUNK;
    outstream1.write(reinterpret_cast<const char*>(recs.data() + written),
                     std::streamsize(block) * std::streamsize(sizeof(mpcdet)));
    if(!outstream1) {
      cerr << "ERROR: write failed after " << written << " records\n";
      return(1);
    }
    written += block;
  }
  outstream1.close();
  cout << "Wrote " << written << " records (" << double(written)*sizeof(mpcdet)/1.0e9
       << " GB) to " << outfile << "\n";

  // Write the text sidecar header for validation/provenance.
  string hdrfile = outfile + ".hdr";
  ofstream hdrstream(hdrfile);
  if(!hdrstream) {
    cerr << "WARNING: can't write header sidecar " << hdrfile << "\n";
  } else {
    hdrstream << "mpcat_index v1\n";
    hdrstream << "nrec " << nrec << "\n";
    hdrstream << "recsize " << sizeof(mpcdet) << "\n";
    hdrstream << fixed << setprecision(8);
    hdrstream << "tmin " << tmin << "\n";
    hdrstream << "tmax " << tmax << "\n";
    hdrstream << "kept " << total_kept << "\n";
    hdrstream << "skipped " << total_skipped << "\n";
    for(i = 0; i < long(source_files.size()); i++) hdrstream << "source " << source_files[i] << "\n";
    hdrstream.close();
    cout << "Wrote header sidecar " << hdrfile << "\n";
  }

  cout << "mpcat_index: done.\n";
  return(0);
}
