//
// Created by delus on 3/29/24.
//

#ifndef TDSPP_OUTPUT_H
#define TDSPP_OUTPUT_H

#include <string>

using namespace std;

struct Output {
    bool md = 0; // Flag for MD state, default is false
    bool ddd = 0; // Flag for DDD state, default is false
    int bpExplored = 0; // Number of BP explored
    double runtime = 0; // Total runtime
    double add_runtime = 0; // Additional runtime
    double sp_runtime = 0; // SP runtime
    int arc_total = 0; // Total number of arcs in solution
    int subPathTotal = 1; // Total number of subpaths in solution
    double optVal = 0; // Optimal objective value
    int iter = 0; // Number of iterations

    void writeOutputCSV(string &filename);
};

#endif //TDSPP_OUTPUT_H
