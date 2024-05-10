//
// Created by delus on 3/29/24.
//

#include "Output.h"

#include <fstream>

void Output::writeOutputCSV(string &filename) {
    ofstream myFile(filename, ofstream::app);
    myFile << this->md << ',';
    myFile << this->ddd << ',';
    myFile << this->bpExplored << ',';
    myFile << this->runtime << ',';
    myFile << this->add_runtime << ',';
    myFile << this->sp_runtime << ',';
    myFile << this->arc_total << ',';
    myFile << this->subPathTotal << ',';
    myFile << this->iter << ',';
}