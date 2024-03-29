//
// Created by delus on 3/29/24.
//

#include "SummaryOutput.h"

void SummaryOutput::writeSOutputCSV() {
    ofstream myFile(filename, ofstream::app);
    myFile << n << ',';
    myFile << T << ',';
    myFile << gtype << ',';
    myFile << ttype << ',';
    myFile << md << ',';
    myFile << avg_bpExplored << ',';
    myFile << bp_total << ',';
    myFile << percent_bp << ',';
    myFile << pow(gavg_percent_time, 1.0 / numSeed) << ',';
    myFile << avg_add_runtime << ',';
    myFile << avg_sp_runtime << ',';
    myFile << avg_arc << ',';
    myFile << avg_subpath << ',';
    myFile << avg_runtime << ',';
    myFile << avg_enumtime << ',';
    myFile << percent_time << ',';
    myFile << iter << endl;
    myFile.close();
    return;
}

void SummaryOutput::updateSummaryOutput(Output &output) {
    if (output.ddd) {
        bpExplored += output.bpExplored;
        runtime += output.runtime;
        total_add_runtime += output.add_runtime;
        total_sp_runtime += output.sp_runtime;
        arc_total += output.arc_total;
        subpath_total += output.subPathTotal;
        gavg_percent_time *= output.runtime;
        iter += output.iter;
    } else {
        bp_total = output.bpExplored;
        enumtime += output.runtime;
        gavg_percent_time /= output.runtime;
    }
    return;
}
void SummaryOutput::calcSummaryOutput() {
    avg_arc = (double) arc_total / numSeed;
    avg_bpExplored = (double) bpExplored / numSeed;
    avg_enumtime = enumtime / numSeed;
    avg_runtime = runtime / numSeed;
    avg_add_runtime = total_add_runtime / numSeed;
    avg_sp_runtime = total_sp_runtime / numSeed;
    avg_subpath = (double) subpath_total / numSeed;
    percent_time = 100 * avg_runtime / avg_enumtime;
    percent_bp = 100 * (double) avg_bpExplored / bp_total;
    return;
}

SummaryOutput::SummaryOutput(string filename, int n, int T, int gtype, int ttype, int numseed, bool md) {
    this->filename = filename;
    this->n = n;
    this->T = T;
    this->gtype = gtype;
    this->ttype = ttype;
    this->numSeed = numseed;
    this->md = md;
}