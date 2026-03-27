// make_hypmask.cpp: March 2026
// Given a heliolinc cluster summary file (coarse run) and a fine hypothesis file,
// output the indices of fine hypotheses that are "near" any found cluster in
// hypothesis space. Used for coarse-to-fine adaptive hypothesis pruning.
//
// Usage:
//   make_hypmask -clust <clustfile> -hyp <hypfile> [-window <float>] -out <outfile>
//
// -clust <file>:   cluster summary CSV from a coarse heliolinc run
// -hyp <file>:     fine hypothesis file (hlradhyp format: r_AU r_dot r_ddot per line)
// -window <float>: fractional search radius in hypothesis space (default 0.5)
//                  a fine hypothesis i is active if there exists a cluster with
//                  |r_fine - r_clust|/r_clust < window  AND
//                  |rdot_fine - rdot_clust| / max(|rdot_clust|, 0.001) < window
// -out <file>:     output file (sorted active indices, one per line)

#include "solarsyst_dyn_geo01.h"
#include <set>
#include <cmath>

static void show_usage()
{
  cerr << "Usage: make_hypmask -clust cluster_summary_file -hyp fine_hypothesis_file\n";
  cerr << "         [-window search_radius_fraction] -out output_index_file\n";
  cerr << "\n";
  cerr << "  -clust <file>:   cluster summary CSV from a coarse heliolinc run\n";
  cerr << "  -hyp <file>:     fine hypothesis file (same format as heliodist input)\n";
  cerr << "  -window <float>: fractional search radius (default 0.5)\n";
  cerr << "                   a fine hypothesis is active if there exists a cluster with\n";
  cerr << "                     |r_fine - r_clust|/r_clust < window  AND\n";
  cerr << "                     |rdot_fine - rdot_clust|/max(|rdot_clust|,0.001) < window\n";
  cerr << "  -out <file>:     output file (sorted active indices, one per line)\n";
}

int main(int argc, char *argv[])
{
  string clustfile, hypfile, outfile;
  double window = 0.5;
  int i = 1;

  if(argc < 7) {
    cerr << "Too few arguments.\n";
    show_usage();
    return(1);
  }

  while(i < argc) {
    if(string(argv[i]) == "-clust" || string(argv[i]) == "--clust") {
      if(i+1 < argc) {
        clustfile = argv[++i];
        i++;
      } else {
        cerr << "ERROR: -clust requires a filename argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-hyp" || string(argv[i]) == "--hyp") {
      if(i+1 < argc) {
        hypfile = argv[++i];
        i++;
      } else {
        cerr << "ERROR: -hyp requires a filename argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-window" || string(argv[i]) == "--window") {
      if(i+1 < argc) {
        window = stod(argv[++i]);
        i++;
      } else {
        cerr << "ERROR: -window requires a numeric argument\n";
        show_usage();
        return(1);
      }
    } else if(string(argv[i]) == "-out" || string(argv[i]) == "--out") {
      if(i+1 < argc) {
        outfile = argv[++i];
        i++;
      } else {
        cerr << "ERROR: -out requires a filename argument\n";
        show_usage();
        return(1);
      }
    } else {
      cerr << "Warning: unrecognized keyword or argument " << argv[i] << "\n";
      i++;
    }
  }

  if(clustfile.empty()) {
    cerr << "ERROR: cluster summary file (-clust) is required\n";
    show_usage();
    return(1);
  }
  if(hypfile.empty()) {
    cerr << "ERROR: fine hypothesis file (-hyp) is required\n";
    show_usage();
    return(1);
  }
  if(outfile.empty()) {
    cerr << "ERROR: output file (-out) is required\n";
    show_usage();
    return(1);
  }

  // Read cluster summary file
  vector<hlclust> clustvec;
  int status = read_clustersum_file(clustfile, clustvec, 0);
  if(status != 0) {
    cerr << "ERROR: read_clustersum_file failed with status " << status << " on file " << clustfile << "\n";
    return(1);
  }
  cout << "Read " << clustvec.size() << " clusters from " << clustfile << "\n";

  if(clustvec.empty()) {
    cerr << "WARNING: no clusters found in " << clustfile << " -- output will be empty\n";
  }

  // Read fine hypothesis file
  vector<hlradhyp> radhyp;
  status = read_radhyp_file(hypfile, radhyp, 0);
  if(status != 0) {
    cerr << "ERROR: read_radhyp_file failed with status " << status << " on file " << hypfile << "\n";
    return(1);
  }
  long nhyp = radhyp.size();
  cout << "Read " << nhyp << " fine hypotheses from " << hypfile << "\n";

  if(nhyp <= 0) {
    cerr << "ERROR: fine hypothesis file is empty\n";
    return(1);
  }

  // For each cluster, find all fine hypotheses within the window
  set<long> active_set;

  for(long ci = 0; ci < long(clustvec.size()); ci++) {
    double r_clust    = clustvec[ci].heliohyp0; // AU
    double rdot_clust = clustvec[ci].heliohyp1; // AU/day

    double rdot_scale = max(fabs(rdot_clust), 0.001);

    for(long hi = 0; hi < nhyp; hi++) {
      double r_fine    = radhyp[hi].HelioRad;
      double rdot_fine = radhyp[hi].R_dot;

      double dr   = fabs(r_fine - r_clust) / r_clust;
      double drdot = fabs(rdot_fine - rdot_clust) / rdot_scale;

      if(dr < window && drdot < window) {
        active_set.insert(hi);
      }
    }
  }

  long nactive = active_set.size();
  cout << "Active hypotheses: " << nactive << " out of " << nhyp
       << " (" << fixed << setprecision(1) << (100.0 * nactive / nhyp) << "%)\n";

  // Write output file
  ofstream outstream(outfile);
  if(!outstream.is_open()) {
    cerr << "ERROR: cannot open output file: " << outfile << "\n";
    return(1);
  }
  for(long idx : active_set) {
    outstream << idx << "\n";
  }
  outstream.close();
  cout << "Wrote " << nactive << " active hypothesis indices to " << outfile << "\n";

  return(0);
}
