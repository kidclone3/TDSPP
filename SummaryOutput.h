//
// Created by delus on 3/29/24.
//

#ifndef TDSPP_SUMMARYOUTPUT_H
#define TDSPP_SUMMARYOUTPUT_H

#include <string>
#include <fstream>
#include <math.h>
#include "Output.h"
using namespace std;

class SummaryOutput {
public:
    bool md;
    int n;
    int T;
    int gtype;
    int ttype;
    int numSeed;
    int bpExplored = 0;
    double runtime = 0;
    double enumtime = 0;
    int arc_total = 0;
    int subpath_total = 0;
    double avg_bpExplored = 0;
    int bp_total = 0;
    double percent_bp = 0;
    double avg_runtime = 0;
    double avg_enumtime = 0;
    double percent_time = 0;
    double gavg_percent_time = 1;
    double avg_arc = 0;
    double avg_subpath = 0;
    double total_add_runtime = 0;
    double total_sp_runtime = 0;
    double avg_add_runtime = 0;
    double avg_sp_runtime = 0;
    int iter = 0;
    string filename;

    SummaryOutput(string filename, int n, int T, int gtype, int ttype, int numseed, bool md);

    void updateSummaryOutput(Output &output);

    void calcSummaryOutput();

    void writeSOutputCSV();
};



#endif //TDSPP_SUMMARYOUTPUT_H
