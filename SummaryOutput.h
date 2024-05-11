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
    bool md; // Flag for MD state, default is false
    int n; // Number of nodes
    int T; // time horizon
    int gtype; // Graph type
    int ttype; // Travel time type
    int numSeed; // seed number
    int bpExplored = 0; // Number of BP explored
    double runtime = 0; // Total runtime for DDD
    double enumtime = 0; // Total runtime for enumerate method
    int arc_total = 0; // Total number of arcs in solution
    int subpath_total = 0; // Total number of subpaths in solution
    double avg_bpExplored = 0; // Average number of BP explored
    int bp_total = 0; // Total number of BP explored
    double percent_bp = 0; // Percentage of BP explored
    double avg_runtime = 0; // Average runtime
    double avg_enumtime = 0; // Average runtime for enumerate method
    double percent_time = 0; // Percentage of runtime vs enumerate method
    double gavg_percent_time = 1; // Geometric mean of percentage of runtime vs enumerate method
    double avg_arc = 0; // Average number of arcs in solution
    double avg_subpath = 0; // Average number of subpaths in solution
    double total_add_runtime = 0; // Total additional runtime
    double total_sp_runtime = 0; // Total SP runtime
    double avg_add_runtime = 0; // Average additional runtime
    double avg_sp_runtime = 0; // Average SP runtime
    int iter = 0; // Number of iterations
    string filename;

    SummaryOutput(string filename, int n, int T, int gtype, int ttype, int numseed, bool md);

    void updateSummaryOutput(Output &output);

    void calcSummaryOutput();

    void writeSOutputCSV();
};



#endif //TDSPP_SUMMARYOUTPUT_H
