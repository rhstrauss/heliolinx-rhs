// 2026-04-28: dedup_across_windows: Cross-window idstring-based linkage deduplication
// for the heliolinx suite.
//
// Reads per-window (pairdets, purify_sum, purify_c2d) tuples listed in an inlist file,
// builds the global cluster pool keyed on ATLAS idstrings, sorts by quality metric,
// and greedily accepts clusters whose idstrings do not overlap with already-accepted ones.
//
// Outputs:
//   -out_csv   : combined parseout.csv  (same format as parse_clust2det)
//   -out_mpc80 : combined parseout.mpc80 (same format as parse_clust2det_MPC80)
//   -out_summary (optional): one line per accepted cluster with dedup provenance

#include "solarsyst_dyn_geo01.h"
#include "cmath"
#include <unordered_set>

static void show_usage()
{
  cerr << "Usage: dedup_across_windows -inlist <file> -out_csv <file> -out_mpc80 <file>\n";
  cerr << "                            [-out_summary <file>] [-trackdiv <hours>]\n";
  cerr << "                            [-min_shared_idstrings <int>]\n";
  cerr << "\n";
  cerr << "inlist format: one line per window, 3 whitespace-separated paths:\n";
  cerr << "  <pairdets_path> <purify_sum_path> <purify_c2d_path>\n";
}

// Per-window cluster record: holds everything needed to write outputs and do dedup
struct winclust {
  int window_id;            // 0-based index into inlist
  long local_clusternum;   // clusternum as it appears in that window's purify_sum
  hlclust cinfo;           // full cluster summary record
  vector<hldet> dets;      // time-sorted detection records (from pairdets)
  vector<string> idstrings; // idstrings of all detections (parallel to dets)
};

// Comparator: sort by metric descending, then astromRMS ascending on ties
struct better_winclust {
  inline bool operator()(const winclust &a, const winclust &b) const {
    if(a.cinfo.metric != b.cinfo.metric) return a.cinfo.metric > b.cinfo.metric;
    return a.cinfo.astromRMS < b.cinfo.astromRMS;
  }
};

// Write the parse_clust2det.cpp-style header+detection block for one cluster
static void write_csv_cluster(ofstream &out, const hlclust &ci, const vector<hldet> &clustvec,
                               long output_clusternum, double nightstep)
{
  long i;
  double angvel, crosstrack, alongtrack, PA, poleRA, poleDec;
  angvel = crosstrack = alongtrack = PA = poleRA = poleDec = 0.0;
  double arc, timespan;
  arc = timespan = 0.0;
  double avg_det_qual = 0.0;
  long max_known_obj = 0;
  double magrange, magmean, magrms, minGCR, maxGCR;
  magrange = magmean = magrms = minGCR = maxGCR = 0.0;
  double crossquad, alongquad, GCRA, GCDec, meantime;
  crossquad = alongquad = GCRA = GCDec = meantime = 0.0;
  double min_nightstep, max_nightstep;
  min_nightstep = max_nightstep = 0.0;

  vector<double> angvelvec, GCRvec, PAvec, arcvec, timespanvec, nightstepvec, magvec;
  vector<double> crossvec, alongvec, timevec, quadfitvec;
  vector<hldet> trackvec;

  if(clustvec.empty()) return;

  for(i = 0; i < long(clustvec.size()); i++) {
    avg_det_qual += double(clustvec[i].det_qual);
    if(clustvec[i].known_obj > max_known_obj) max_known_obj = clustvec[i].known_obj;
    if(clustvec[i].mag > 0.0) magvec.push_back(clustvec[i].mag);
  }
  avg_det_qual /= double(clustvec.size());

  greatcircfit(clustvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
  for(i = 0; i < long(clustvec.size()); i++) meantime += clustvec[i].MJD / double(clustvec.size());
  crossvec = alongvec = timevec = {};
  for(i = 0; i < long(clustvec.size()); i++) {
    poleswitch02(clustvec[i].RA, clustvec[i].Dec, poleRA, poleDec, 90.0, GCRA, GCDec);
    alongvec.push_back(GCRA);
    crossvec.push_back(GCDec);
    timevec.push_back(clustvec[i].MJD - meantime);
  }
  double unwrapping = 0.0;
  for(i = 1; i < long(clustvec.size()); i++) {
    if(alongvec[i] > alongvec[i-1] + 180.0) unwrapping = -360.0;
    else if(alongvec[i] < alongvec[i-1] - 180.0) unwrapping = +360.0;
    alongvec[i] += unwrapping;
  }
  quadfitvec = {};
  polyfit01(alongvec, timevec, int(clustvec.size()), 2, quadfitvec);
  alongquad = quadfitvec[2];
  quadfitvec = {};
  polyfit01(crossvec, timevec, int(clustvec.size()), 2, quadfitvec);
  crossquad = quadfitvec[2];

  // Tracklet analysis (mirrors parse_clust2det exactly)
  trackvec = {};
  nightstepvec = {};
  trackvec.push_back(clustvec[0]);
  angvelvec = GCRvec = PAvec = timespanvec = arcvec = {};
  for(i = 1; i < long(clustvec.size()); i++) {
    if((clustvec[i].MJD - clustvec[i-1].MJD) >= nightstep) {
      nightstepvec.push_back(clustvec[i].MJD - clustvec[i-1].MJD);
    }
    if((clustvec[i].MJD - clustvec[i-1].MJD) < nightstep) {
      trackvec.push_back(clustvec[i]);
    } else {
      long tracknum = trackvec.size();
      if(tracknum > 1) {
        greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
        angvelvec.push_back(angvel);
        if(tracknum > 2) GCRvec.push_back(sqrt(DSQUARE(crosstrack) + DSQUARE(alongtrack)));
        PAvec.push_back(PA);
        timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
        arc = timespan * angvel;
        arcvec.push_back(arc * 3600.0);
        timespanvec.push_back(timespan * 24.0);
        trackvec = {};
        trackvec.push_back(clustvec[i]);
      } else {
        angvelvec.push_back(-1.0);
        PAvec.push_back(-999.0);
        arcvec.push_back(0.0);
        timespanvec.push_back(0.0);
        trackvec = {};
        trackvec.push_back(clustvec[i]);
      }
    }
  }
  {
    long tracknum = trackvec.size();
    if(tracknum > 1) {
      greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
      angvelvec.push_back(angvel);
      if(tracknum > 2) GCRvec.push_back(sqrt(DSQUARE(crosstrack) + DSQUARE(alongtrack)));
      PAvec.push_back(PA);
      timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
      arc = timespan * angvel;
      arcvec.push_back(arc * 3600.0);
      timespanvec.push_back(timespan * 24.0);
      trackvec = {};
    } else {
      angvelvec.push_back(-1.0);
      PAvec.push_back(-999.0);
      arcvec.push_back(0.0);
      timespanvec.push_back(0.0);
      trackvec = {};
    }
  }
  long tracknum = angvelvec.size();
  sort(angvelvec.begin(), angvelvec.end());
  sort(PAvec.begin(), PAvec.end());
  sort(timespanvec.begin(), timespanvec.end());
  sort(arcvec.begin(), arcvec.end());

  if(long(GCRvec.size()) < 1) minGCR = maxGCR = 0.0;
  else if(long(GCRvec.size()) == 1) minGCR = maxGCR = GCRvec[0];
  else {
    sort(GCRvec.begin(), GCRvec.end());
    minGCR = GCRvec[0];
    maxGCR = GCRvec[GCRvec.size()-1];
  }

  if(!nightstepvec.empty()) {
    sort(nightstepvec.begin(), nightstepvec.end());
    min_nightstep = nightstepvec[0];
    max_nightstep = nightstepvec[nightstepvec.size()-1];
  }

  if(magvec.empty()) {
    magmean = 0.0; magrms = magrange = 99.9;
  } else if(magvec.size() == 1) {
    magmean = magvec[0]; magrms = magrange = 99.9;
  } else if(magvec.size() <= 5) {
    sort(magvec.begin(), magvec.end());
    dmeanrms01(magvec, &magmean, &magrms);
    magrange = magvec[magvec.size()-1] - magvec[0];
  } else {
    sort(magvec.begin(), magvec.end());
    dmeanrms01(magvec, &magmean, &magrms);
    magrange = magvec[magvec.size()-2] - magvec[1];
  }

  out << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_incl,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,stringID,min_nightstep,max_nightstep,magmean,magrms,magrange,rating,crossaccel,alongaccel,totalaccel\n";
  out << fixed << setprecision(3) << output_clusternum << "," << ci.posRMS << "," << ci.velRMS << "," << ci.totRMS << ",";
  out << fixed << setprecision(4) << ci.astromRMS << ",";
  out << fixed << setprecision(6) << ci.timespan << "," << ci.uniquepoints << "," << ci.obsnights << "," << ci.metric << ",";
  out << fixed << setprecision(6) << ci.orbit_a << "," << ci.orbit_e << "," << ci.orbit_incl << "," << ci.orbit_MJD << ",";
  out << fixed << setprecision(1) << ci.orbitX << "," << ci.orbitY << "," << ci.orbitZ << ",";
  out << fixed << setprecision(4) << ci.orbitVX << "," << ci.orbitVY << "," << ci.orbitVZ << "," << ci.orbit_eval_count << ",";
  out << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
  out << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << minGCR << "," << maxGCR << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << clustvec[0].idstring << "," << min_nightstep << "," << max_nightstep << "," << magmean << "," << magrms << "," << magrange << "," << ci.rating << "," << crossquad << "," << alongquad << "," << sqrt(crossquad*crossquad + alongquad*alongquad) << "\n";

  out << "#MJD,RA,Dec,mag,trail_len,trail_PA,sigmag,sig_across,sig_along,image,idstring,band,obscode,known_obj,det_qual,clusternum\n";
  for(i = 0; i < long(clustvec.size()); i++) {
    out << fixed << setprecision(7) << clustvec[i].MJD << "," << clustvec[i].RA << "," << clustvec[i].Dec << ",";
    out << fixed << setprecision(4) << clustvec[i].mag << ",";
    out << fixed << setprecision(2) << clustvec[i].trail_len << "," << clustvec[i].trail_PA << ",";
    out << fixed << setprecision(4) << clustvec[i].sigmag << ",";
    out << fixed << setprecision(3) << clustvec[i].sig_across << "," << clustvec[i].sig_along << ",";
    out << clustvec[i].image << "," << clustvec[i].idstring << "," << clustvec[i].band << ",";
    out << clustvec[i].obscode << "," << clustvec[i].known_obj << ",";
    out << clustvec[i].det_qual << "," << output_clusternum << "\n";
  }
}

// Write the parse_clust2det_MPC80.cpp-style block for one cluster
static void write_mpc80_cluster(ofstream &out, const hlclust &ci, const vector<hldet> &clustvec,
                                 long output_clusternum, double nightstep)
{
  long i;
  double angvel, crosstrack, alongtrack, PA, poleRA, poleDec;
  angvel = crosstrack = alongtrack = PA = poleRA = poleDec = 0.0;
  double arc, timespan;
  arc = timespan = 0.0;
  double avg_det_qual = 0.0;
  long max_known_obj = 0;
  double magrange, magmean, magrms, minGCR, maxGCR;
  magrange = magmean = magrms = minGCR = maxGCR = 0.0;
  double crossquad, alongquad, GCRA, GCDec, meantime;
  crossquad = alongquad = GCRA = GCDec = meantime = 0.0;
  double min_nightstep, max_nightstep;
  min_nightstep = max_nightstep = 0.0;

  vector<double> angvelvec, GCRvec, PAvec, arcvec, timespanvec, nightstepvec, magvec;
  vector<double> crossvec, alongvec, timevec, quadfitvec;
  vector<hldet> trackvec;

  if(clustvec.empty()) return;

  for(i = 0; i < long(clustvec.size()); i++) {
    avg_det_qual += double(clustvec[i].det_qual);
    if(clustvec[i].known_obj > max_known_obj) max_known_obj = clustvec[i].known_obj;
    if(clustvec[i].mag > 0.0) magvec.push_back(clustvec[i].mag);
  }
  avg_det_qual /= double(clustvec.size());

  greatcircfit(clustvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
  for(i = 0; i < long(clustvec.size()); i++) meantime += clustvec[i].MJD / double(clustvec.size());
  crossvec = alongvec = timevec = {};
  for(i = 0; i < long(clustvec.size()); i++) {
    poleswitch02(clustvec[i].RA, clustvec[i].Dec, poleRA, poleDec, 90.0, GCRA, GCDec);
    alongvec.push_back(GCRA);
    crossvec.push_back(GCDec);
    timevec.push_back(clustvec[i].MJD - meantime);
  }
  double unwrapping = 0.0;
  for(i = 1; i < long(clustvec.size()); i++) {
    if(alongvec[i] > alongvec[i-1] + 180.0) unwrapping = -360.0;
    else if(alongvec[i] < alongvec[i-1] - 180.0) unwrapping = +360.0;
    alongvec[i] += unwrapping;
  }
  quadfitvec = {};
  polyfit01(alongvec, timevec, int(clustvec.size()), 2, quadfitvec);
  alongquad = quadfitvec[2];
  quadfitvec = {};
  polyfit01(crossvec, timevec, int(clustvec.size()), 2, quadfitvec);
  crossquad = quadfitvec[2];

  trackvec = {};
  nightstepvec = {};
  trackvec.push_back(clustvec[0]);
  angvelvec = GCRvec = PAvec = timespanvec = arcvec = {};
  for(i = 1; i < long(clustvec.size()); i++) {
    if((clustvec[i].MJD - clustvec[i-1].MJD) >= nightstep) {
      nightstepvec.push_back(clustvec[i].MJD - clustvec[i-1].MJD);
    }
    if((clustvec[i].MJD - clustvec[i-1].MJD) < nightstep) {
      trackvec.push_back(clustvec[i]);
    } else {
      long tracknum = trackvec.size();
      if(tracknum > 1) {
        greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
        angvelvec.push_back(angvel);
        if(tracknum > 2) GCRvec.push_back(sqrt(DSQUARE(crosstrack) + DSQUARE(alongtrack)));
        PAvec.push_back(PA);
        timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
        arc = timespan * angvel;
        arcvec.push_back(arc * 3600.0);
        timespanvec.push_back(timespan * 24.0);
        trackvec = {};
        trackvec.push_back(clustvec[i]);
      } else {
        angvelvec.push_back(-1.0);
        PAvec.push_back(-999.0);
        arcvec.push_back(0.0);
        timespanvec.push_back(0.0);
        trackvec = {};
        trackvec.push_back(clustvec[i]);
      }
    }
  }
  {
    long tracknum = trackvec.size();
    if(tracknum > 1) {
      greatcircfit(trackvec, poleRA, poleDec, angvel, PA, crosstrack, alongtrack);
      angvelvec.push_back(angvel);
      if(tracknum > 2) GCRvec.push_back(sqrt(DSQUARE(crosstrack) + DSQUARE(alongtrack)));
      PAvec.push_back(PA);
      timespan = trackvec[trackvec.size()-1].MJD - trackvec[0].MJD;
      arc = timespan * angvel;
      arcvec.push_back(arc * 3600.0);
      timespanvec.push_back(timespan * 24.0);
      trackvec = {};
    } else {
      angvelvec.push_back(-1.0);
      PAvec.push_back(-999.0);
      arcvec.push_back(0.0);
      timespanvec.push_back(0.0);
      trackvec = {};
    }
  }
  long tracknum = angvelvec.size();
  sort(angvelvec.begin(), angvelvec.end());
  sort(PAvec.begin(), PAvec.end());
  sort(timespanvec.begin(), timespanvec.end());
  sort(arcvec.begin(), arcvec.end());

  if(long(GCRvec.size()) < 1) minGCR = maxGCR = 0.0;
  else if(long(GCRvec.size()) == 1) minGCR = maxGCR = GCRvec[0];
  else {
    sort(GCRvec.begin(), GCRvec.end());
    minGCR = GCRvec[0];
    maxGCR = GCRvec[GCRvec.size()-1];
  }

  if(!nightstepvec.empty()) {
    sort(nightstepvec.begin(), nightstepvec.end());
    min_nightstep = nightstepvec[0];
    max_nightstep = nightstepvec[nightstepvec.size()-1];
  }

  if(magvec.empty()) {
    magmean = 0.0; magrms = magrange = 99.9;
  } else if(magvec.size() == 1) {
    magmean = magvec[0]; magrms = magrange = 99.9;
  } else if(magvec.size() <= 5) {
    sort(magvec.begin(), magvec.end());
    dmeanrms01(magvec, &magmean, &magrms);
    magrange = magvec[magvec.size()-1] - magvec[0];
  } else {
    sort(magvec.begin(), magvec.end());
    dmeanrms01(magvec, &magmean, &magrms);
    magrange = magvec[magvec.size()-2] - magvec[1];
  }

  out << "\n#clusternum,posRMS,velRMS,totRMS,astromRMS,timespan,uniquepoints,obsnights,metric,orbit_a,orbit_e,orbit_incl,orbit_MJD,orbitX,orbitY,orbitZ,orbitVX,orbitVY,orbitVZ,orbit_eval_count,avg_det_qual,max_known_obj,minvel,maxvel,minGCR,maxGCR,minpa,maxpa,mintimespan,maxtimespan,minarc,maxarc,stringID,min_nightstep,max_nightstep,magmean,magrms,magrange,rating,crossaccel,alongaccel,totalaccel\n";
  out << fixed << setprecision(3) << output_clusternum << "," << ci.posRMS << "," << ci.velRMS << "," << ci.totRMS << ",";
  out << fixed << setprecision(4) << ci.astromRMS << ",";
  out << fixed << setprecision(6) << ci.timespan << "," << ci.uniquepoints << "," << ci.obsnights << "," << ci.metric << ",";
  out << fixed << setprecision(6) << ci.orbit_a << "," << ci.orbit_e << "," << ci.orbit_incl << "," << ci.orbit_MJD << ",";
  out << fixed << setprecision(1) << ci.orbitX << "," << ci.orbitY << "," << ci.orbitZ << ",";
  out << fixed << setprecision(4) << ci.orbitVX << "," << ci.orbitVY << "," << ci.orbitVZ << "," << ci.orbit_eval_count << ",";
  out << fixed << setprecision(1) << avg_det_qual << "," << max_known_obj << ",";
  out << fixed << setprecision(6) << angvelvec[0] << "," << angvelvec[tracknum-1] << "," << minGCR << "," << maxGCR << "," << PAvec[0] << "," << PAvec[tracknum-1] << "," << timespanvec[0] << "," << timespanvec[tracknum-1] << "," << arcvec[0] << "," << arcvec[tracknum-1] << "," << clustvec[0].idstring << "," << min_nightstep << "," << max_nightstep << "," << magmean << "," << magrms << "," << magrange << "," << ci.rating << "," << crossquad << "," << alongquad << "," << sqrt(crossquad*crossquad + alongquad*alongquad) << "\n";

  out << "#MPC 80-column formatted observations:\n";
  int year, month;
  double day;
  int rahr, ramin;
  int decdeg, decmin;
  double decsec, rasec, Dec;
  string stest1, signstring;
  long j;
  int bandlen;

  for(i = 0; i < long(clustvec.size()); i++) {
    // Temporary ID: 'a' + 6-digit output_clusternum
    if(output_clusternum <= 999999) stest1 = to_string(output_clusternum);
    else stest1 = to_string(output_clusternum % 1000000);
    while(stest1.size() < 6) stest1 = "0" + stest1;
    out << "     a" << stest1 << "  C";
    mjd2mpcdate(clustvec[i].MJD, year, month, day);
    out << year << " ";
    if(month < 10) out << "0";
    out << month << " ";
    day = round(day * 1000000.0) / 1000000.0;
    if(day < 10.0) out << fixed << setprecision(6) << "0";
    out << fixed << setprecision(6) << day;
    // RA
    rahr  = int(clustvec[i].RA / 15.0);
    ramin = int(clustvec[i].RA * 4.0 - double(rahr) * 60.0);
    rasec = clustvec[i].RA * 240.0 - double(rahr) * 3600.0 - double(ramin) * 60.0;
    rasec = round(rasec * 1000.0) / 1000.0;
    // Dec
    if(clustvec[i].Dec >= 0) { signstring = "+"; Dec = clustvec[i].Dec; }
    else                      { signstring = "-"; Dec = -clustvec[i].Dec; }
    decdeg = int(Dec);
    decmin = int(Dec * 60.0 - double(decdeg) * 60.0);
    decsec = Dec * 3600.0 - double(decdeg) * 3600.0 - double(decmin) * 60.0;
    decsec = round(decsec * 100.0) / 100.0;
    if(rahr < 10) out << "0";
    out << rahr << " ";
    if(ramin < 10) out << "0";
    out << ramin << " ";
    if(rasec < 10.0) out << "0";
    out << fixed << setprecision(3) << rasec << signstring;
    if(decdeg < 10) out << "0";
    out << decdeg << " ";
    if(decmin < 10) out << "0";
    out << decmin << " ";
    if(decsec < 10.0) out << "0";
    out << fixed << setprecision(2) << decsec << "         ";
    out << fixed << setprecision(1) << clustvec[i].mag << " " << clustvec[i].band;
    bandlen = j = 0;
    while(j < MINSTRINGLEN && clustvec[i].band[j] != '\0') { bandlen++; j++; }
    for(j = 0; j < 7 - bandlen; j++) out << " ";
    out << clustvec[i].obscode << "\n";
  }
}


int main(int argc, char *argv[])
{
  string inlistfile, out_csv_file, out_mpc80_file, out_summary_file;
  double nightstep = NIGHTSTEP;
  int min_shared_idstrings = 3;
  int verbose = 0;
  int status = 0;
  long i;

  if(argc < 7) {
    show_usage();
    return 1;
  }

  i = 1;
  while(i < argc) {
    if(string(argv[i]) == "-inlist") {
      if(i + 1 < argc) { inlistfile = argv[++i]; i++; }
      else { cerr << "ERROR: -inlist requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-out_csv") {
      if(i + 1 < argc) { out_csv_file = argv[++i]; i++; }
      else { cerr << "ERROR: -out_csv requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-out_mpc80") {
      if(i + 1 < argc) { out_mpc80_file = argv[++i]; i++; }
      else { cerr << "ERROR: -out_mpc80 requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-out_summary") {
      if(i + 1 < argc) { out_summary_file = argv[++i]; i++; }
      else { cerr << "ERROR: -out_summary requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-trackdiv" || string(argv[i]) == "-nightstep") {
      if(i + 1 < argc) { nightstep = stod(argv[++i]) / 24.0; i++; }
      else { cerr << "ERROR: -trackdiv requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-min_shared_idstrings") {
      if(i + 1 < argc) { min_shared_idstrings = stoi(argv[++i]); i++; }
      else { cerr << "ERROR: -min_shared_idstrings requires an argument\n"; show_usage(); return 1; }
    } else if(string(argv[i]) == "-verbose" || string(argv[i]) == "-verb") {
      if(i + 1 < argc) { verbose = stoi(argv[++i]); i++; }
      else { cerr << "ERROR: -verbose requires an argument\n"; show_usage(); return 1; }
    } else {
      cerr << "Warning: unrecognized argument " << argv[i] << "\n";
      i++;
    }
  }

  if(inlistfile.empty()) { cerr << "ERROR: -inlist is required\n"; show_usage(); return 1; }
  if(out_csv_file.empty()) { cerr << "ERROR: -out_csv is required\n"; show_usage(); return 1; }
  if(out_mpc80_file.empty()) { cerr << "ERROR: -out_mpc80 is required\n"; show_usage(); return 1; }

  cout << "dedup_across_windows: reading inlist " << inlistfile << "\n";
  cout << "  out_csv: " << out_csv_file << "\n";
  cout << "  out_mpc80: " << out_mpc80_file << "\n";
  cout << "  min_shared_idstrings: " << min_shared_idstrings << "\n";
  cout << "  nightstep: " << nightstep * 24.0 << " hours\n";

  // ── PHASE 1: read all windows and build global cluster pool ──────────────

  vector<winclust> pool;
  long total_input_clusters = 0;
  int window_id = 0;

  ifstream inlist_stream(inlistfile);
  if(!inlist_stream.is_open()) {
    cerr << "ERROR: could not open inlist file " << inlistfile << "\n";
    return 1;
  }

  string line;
  while(getline(inlist_stream, line)) {
    // Skip blank lines and comments
    if(line.empty()) continue;
    size_t first_nonspace = line.find_first_not_of(" \t\r\n");
    if(first_nonspace == string::npos) continue;
    if(line[first_nonspace] == '#') continue;

    istringstream ss(line);
    string pairdet_path, sum_path, c2d_path;
    ss >> pairdet_path >> sum_path >> c2d_path;
    if(pairdet_path.empty() || sum_path.empty() || c2d_path.empty()) {
      cerr << "WARNING: skipping malformed inlist line: " << line << "\n";
      continue;
    }

    cout << "\nWindow " << window_id << ":\n";
    cout << "  pairdets: " << pairdet_path << "\n";
    cout << "  sum:      " << sum_path << "\n";
    cout << "  c2d:      " << c2d_path << "\n";

    // Read pairdets
    vector<hldet> detvec;
    status = read_pairdet_file(pairdet_path, detvec, verbose);
    if(status != 0) {
      cerr << "ERROR: could not read pairdet file " << pairdet_path << " (status=" << status << ")\n";
      return 1;
    }
    cout << "  read " << detvec.size() << " pairdet entries\n";

    // Build detnum -> idstring map (detnum = index in detvec, 0-based)
    // detvec[i].idstring is directly available
    // (pairdet format has idstring as a field in hldet)

    // Read cluster summary
    vector<hlclust> inclustvec;
    status = read_clustersum_file(sum_path, inclustvec, verbose);
    if(status != 0) {
      cerr << "ERROR: could not read cluster summary file " << sum_path << " (status=" << status << ")\n";
      return 1;
    }
    cout << "  read " << inclustvec.size() << " cluster summary entries\n";

    // Read cluster-to-detection file
    vector<longpair> inclust2det;
    status = read_longpair_file(c2d_path, inclust2det, verbose);
    if(status != 0) {
      cerr << "ERROR: could not read c2d file " << c2d_path << " (status=" << status << ")\n";
      return 1;
    }
    cout << "  read " << inclust2det.size() << " c2d entries\n";

    // Build per-cluster detection lists using parse_clust2det library function
    vector<hldet> cluster_detvec;
    status = parse_clust2det(detvec, inclust2det, cluster_detvec);
    if(status != 0) {
      cerr << "ERROR: parse_clust2det failed (status=" << status << ") for window " << window_id << "\n";
      return 1;
    }

    // Now iterate over clusters in inclustvec, collect detection records
    long detnum = cluster_detvec.size();
    long detct = 0;
    for(long clustct = 0; clustct < long(inclustvec.size()); clustct++) {
      vector<hldet> clustvec;
      if(detct < detnum) {
        if(cluster_detvec[detct].index != clustct) {
          cerr << "ERROR: cluster counting mismatch in window " << window_id
               << " clustct=" << clustct << " detct=" << detct
               << " index=" << cluster_detvec[detct].index << "\n";
          return 2;
        }
        while(detct < detnum && cluster_detvec[detct].index == clustct) {
          clustvec.push_back(cluster_detvec[detct]);
          detct++;
        }
      }
      if(clustvec.empty()) continue;

      // Time-sort
      sort(clustvec.begin(), clustvec.end(), early_hldet());

      winclust wc;
      wc.window_id = window_id;
      wc.local_clusternum = inclustvec[clustct].clusternum;
      wc.cinfo = inclustvec[clustct];
      wc.dets = clustvec;
      for(long p = 0; p < long(clustvec.size()); p++) {
        wc.idstrings.push_back(string(clustvec[p].idstring));
      }
      pool.push_back(wc);
      total_input_clusters++;
    }

    window_id++;
  }
  inlist_stream.close();

  int total_windows = window_id;
  cout << "\nTotal windows loaded: " << total_windows << "\n";
  cout << "Total input clusters: " << total_input_clusters << "\n";

  // ── PHASE 2: sort pool by quality (metric desc, astromRMS asc) ───────────

  sort(pool.begin(), pool.end(), better_winclust());

  // ── PHASE 3: greedy idstring-based dedup ─────────────────────────────────

  unordered_set<string> globally_used;
  vector<long> accepted_pool_indices;   // indices into pool[] for accepted clusters
  vector<long> dup_winner_for;          // accepted index (into accepted_pool_indices)
                                         // that this duplicate maps to, or -1 if accepted
  // We'll track which accepted cluster each duplicate was absorbed by
  // For summary: per accepted cluster, list of (window, local_clusternum) duplicates
  struct accepted_info {
    long pool_idx;
    vector<pair<int,long>> duplicates; // (window_id, local_clusternum)
  };
  vector<accepted_info> accepted;

  for(long k = 0; k < long(pool.size()); k++) {
    const winclust &wc = pool[k];

    // Count overlap with globally used idstrings
    int overlap = 0;
    for(const string &ids : wc.idstrings) {
      if(globally_used.count(ids)) {
        overlap++;
        if(overlap >= min_shared_idstrings) break;
      }
    }

    if(overlap >= min_shared_idstrings) {
      // Duplicate: find which accepted cluster owns the most of our idstrings
      // (We record it under the first accepted cluster that shares any idstring)
      long best_accepted = -1;
      for(long ai = long(accepted.size()) - 1; ai >= 0; ai--) {
        const winclust &awc = pool[accepted[ai].pool_idx];
        int shared = 0;
        unordered_set<string> aset(awc.idstrings.begin(), awc.idstrings.end());
        for(const string &ids : wc.idstrings) {
          if(aset.count(ids)) shared++;
        }
        if(shared >= min_shared_idstrings) {
          best_accepted = ai;
          break;
        }
      }
      // Linear search from back is fine (pool is O(100s) of clusters)
      // If somehow no match found (shouldn't happen), use the last accepted
      if(best_accepted < 0 && !accepted.empty()) best_accepted = long(accepted.size()) - 1;
      if(best_accepted >= 0) {
        accepted[best_accepted].duplicates.push_back({wc.window_id, wc.local_clusternum});
      }
    } else {
      // Accept: mark all idstrings as used
      for(const string &ids : wc.idstrings) {
        globally_used.insert(ids);
      }
      accepted_info ai;
      ai.pool_idx = k;
      accepted.push_back(ai);
    }
  }

  long n_accepted = accepted.size();
  long n_duplicates = total_input_clusters - n_accepted;
  cout << "\nDedup result: " << total_input_clusters << " input -> "
       << n_accepted << " unique, " << n_duplicates << " duplicates collapsed\n";

  // ── PHASE 4: write outputs ────────────────────────────────────────────────

  ofstream csv_out(out_csv_file);
  if(!csv_out.is_open()) {
    cerr << "ERROR: could not open output csv file " << out_csv_file << "\n";
    return 1;
  }
  ofstream mpc80_out(out_mpc80_file);
  if(!mpc80_out.is_open()) {
    cerr << "ERROR: could not open output mpc80 file " << out_mpc80_file << "\n";
    return 1;
  }

  for(long ai = 0; ai < n_accepted; ai++) {
    const winclust &wc = pool[accepted[ai].pool_idx];
    write_csv_cluster(csv_out, wc.cinfo, wc.dets, ai, nightstep);
    write_mpc80_cluster(mpc80_out, wc.cinfo, wc.dets, ai, nightstep);
  }
  csv_out.close();
  mpc80_out.close();
  cout << "Wrote " << n_accepted << " clusters to " << out_csv_file << "\n";
  cout << "Wrote " << n_accepted << " clusters to " << out_mpc80_file << "\n";

  // Optional summary
  if(!out_summary_file.empty()) {
    ofstream sum_out(out_summary_file);
    if(!sum_out.is_open()) {
      cerr << "WARNING: could not open summary file " << out_summary_file << " -- skipping\n";
    } else {
      sum_out << "#output_clusternum,winning_window,winning_local_clusternum,n_detections,metric,astromRMS,obsnights,timespan,n_duplicates_absorbed,duplicate_windows_local_clusternums\n";
      for(long ai = 0; ai < n_accepted; ai++) {
        const winclust &wc = pool[accepted[ai].pool_idx];
        sum_out << ai << "," << wc.window_id << "," << wc.local_clusternum << ","
                << wc.dets.size() << ","
                << fixed << setprecision(6) << wc.cinfo.metric << ","
                << fixed << setprecision(4) << wc.cinfo.astromRMS << ","
                << wc.cinfo.obsnights << ","
                << fixed << setprecision(6) << wc.cinfo.timespan << ","
                << accepted[ai].duplicates.size();
        for(const auto &dup : accepted[ai].duplicates) {
          sum_out << ",win" << dup.first << "_c" << dup.second;
        }
        sum_out << "\n";
      }
      sum_out.close();
      cout << "Wrote summary to " << out_summary_file << "\n";
    }
  }

  // Print some sample dedup groups for diagnostics
  cout << "\n--- Sample duplicate groups (first 5 with duplicates) ---\n";
  int shown = 0;
  for(long ai = 0; ai < n_accepted && shown < 5; ai++) {
    if(accepted[ai].duplicates.empty()) continue;
    const winclust &wc = pool[accepted[ai].pool_idx];
    cout << "  Output clust " << ai << ": winner=win" << wc.window_id
         << "_c" << wc.local_clusternum
         << " metric=" << fixed << setprecision(3) << wc.cinfo.metric
         << " ndets=" << wc.dets.size()
         << " dups:";
    for(const auto &dup : accepted[ai].duplicates) {
      cout << " win" << dup.first << "_c" << dup.second;
    }
    cout << "\n";
    shown++;
  }
  if(shown == 0) cout << "  (no duplicate groups found)\n";

  // Verify: check that no idstring appears more than once across output clusters
  // (expensive check, only in verbose mode)
  if(verbose >= 1) {
    cout << "\nVerbose: verifying no duplicate idstrings in output...\n";
    unordered_map<string,long> ids_to_cluster;
    bool dup_found = false;
    for(long ai = 0; ai < n_accepted; ai++) {
      const winclust &wc = pool[accepted[ai].pool_idx];
      for(const string &ids : wc.idstrings) {
        if(ids_to_cluster.count(ids)) {
          cerr << "  DUPLICATE idstring in output: " << ids
               << " in output clusters " << ids_to_cluster[ids] << " and " << ai << "\n";
          dup_found = true;
        } else {
          ids_to_cluster[ids] = ai;
        }
      }
    }
    if(!dup_found) cout << "  OK: no duplicate idstrings in output.\n";
    else cerr << "  ERROR: duplicate idstrings found in output -- check min_shared_idstrings setting.\n";
  }

  cout << "\ndedup_across_windows complete.\n";
  return 0;
}
