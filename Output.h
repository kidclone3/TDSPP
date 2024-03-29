//
// Created by delus on 3/29/24.
//

#ifndef TDSPP_OUTPUT_H
#define TDSPP_OUTPUT_H
#include <string>
#include <fstream>
using namespace std;

struct Output {
    bool md = 0;
    bool ddd = 0;
    int bpExplored = 0;
    double runtime = 0;
    double add_runtime = 0;
    double sp_runtime = 0;
    int arc_total = 0;
    int subPathTotal = 1;
    double optVal = 0;
    int iter = 0;

    void writeOutputCSV(string &filename);
};

#endif //TDSPP_OUTPUT_H
