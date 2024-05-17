// TDSPP.cpp : This file contains the 'main' function. Program execution begins and ends there.
// TO DO: add more breakpoints associated with minimums (not just THE minimum).


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

// Import self-defined classes:
#include "Output.h"
#include "SummaryOutput.h"
#include "TEN.h"

using namespace std;
typedef pair<int, int> TypeArc;
typedef pair<int, int> TypeBP;

void writeOutputCSVHeader(const string &filename) {
    ofstream myFile(filename, ofstream::app);
    myFile << "Is MD?" << ',';
    myFile << "Is DDD?" << ',';
    myFile << "BP Explored" << ',';
    myFile << "Run-Time" << ',';
    myFile << "Add Run-Time" << ',';
    myFile << "SP Run-Time" << ',';
    myFile << "#Arcs in Soln" << ',';
    myFile << "#Subpaths in Soln" << ',';
    myFile << "Iterations" << ',';
    myFile << "Opt Obj" << endl;
    myFile.close();
    return;
}

void writeSOutputCSVHeader(const string &filename) {
    ofstream myFile(filename, ofstream::app);
    myFile << "n" << ',';
    myFile << "T" << ',';
    myFile << "g-type" << ',';
    myFile << "t-type" << ',';
    myFile << "md" << ',';
    myFile << "Avg BP Explored" << ',';
    myFile << "Total BP" << ',';
    myFile << "BP Explored (%)" << ',';
    myFile << "Geom Run-Time (%)" << ',';
    myFile << "Add Run-Time (ms)" << ',';
    myFile << "SP Run-Time (ms)" << ',';
    myFile << "Avg #Arcs in Soln" << ',';
    myFile << "Avg #Subpaths in Soln" << ',';
    myFile << "Avg DDD Run-Time (ms)" << ',';
    myFile << "Avg Enum Run-Time (ms)" << ',';
    myFile << "Run-Time (%)" << ',';
    myFile << "Iterations" << endl;
    myFile.close();
    return;
}

void runTest(vector<int> ns, vector<int> endTs, vector<int> gTypes, vector<int> tTypes, vector<int> seeds,
             int common_startTime = 0, bool md = 1, bool mtt = 1, bool mdEnum = 1, bool mttEnum = 1) {
    /*
     * ns: vector of number of nodes
     * endTs: vector of end times
     * gTypes: vector of graph types
     * tTypes: vector of travel time types
     * seeds: vector of seeds
     * common_startTime: common start time
     * md: flag for MD: Minimum Duration
     * mtt: flag for MTT: Minimum Travel Time
     * mdEnum: Default value of 1 (true). If true, the function will run the minimum duration enumeration (mdEnum) test.
     * mttEnum: Default value of 1 (true). If true, the function will run the minimum travel time enumeration (mttEnum) test.
     */
//    string md_filename = "Results/experimentsMDMED.csv";
//    string mtt_filename = "Results/MTTBPMMED.csv";

    string md_filename = "Results/experimentsMDMIN.csv";
    string mtt_filename = "Results/MTTBPMMIN.csv";
    if (md) writeSOutputCSVHeader(md_filename);
    if (mtt) writeSOutputCSVHeader(mtt_filename);
    for (int n: ns) {
        for (int common_endTime: endTs) {
            for (int gtype: gTypes) {
                for (int ttype: tTypes) {
                    SummaryOutput myMDOutput(md_filename, n, common_endTime - common_startTime, gtype, ttype,
                                             seeds.size(), 1);
                    SummaryOutput myMTTOutput(mtt_filename, n, common_endTime - common_startTime, gtype, ttype,
                                              seeds.size(),
                                              0);
//                    string filename = "Results/experimentsMDMEDlog.csv";
                    string filename = "Results/experimentsMDMINlog.csv";
                    //string filename = "Results/n" + to_string(n) + "T" + to_string(common_endTime - common_startTime) + "g" + to_string(gtype) + "t" + to_string(ttype) + ".csv";
                    writeOutputCSVHeader(filename);
                    for (int seed: seeds) {
                        if (md) {
                            TEN my_ten(n, common_endTime, gtype, ttype, seed);
                            Output output = my_ten.findMD();
                            output.writeOutputCSV(filename);
                            myMDOutput.updateSummaryOutput(output);
                        }
                        if (mdEnum) {
                            TEN my_ten(n, common_endTime, gtype, ttype, seed);
                            Output output = my_ten.findEnumMD();
                            output.writeOutputCSV(filename);
                            myMDOutput.updateSummaryOutput(output);
                        }
                        if (mtt) {
                            TEN my_ten(n, common_endTime, gtype, ttype, seed);
                            Output output = my_ten.findMTT();
                            output.writeOutputCSV(filename);
                            myMTTOutput.updateSummaryOutput(output);
                        }
                        if (mttEnum) {
                            TEN my_ten(n, common_endTime, gtype, ttype, seed);
                            Output output = my_ten.findEnumMTT();
                            output.writeOutputCSV(filename);
                            myMTTOutput.updateSummaryOutput(output);
                        }
                    }
                    if (md) myMDOutput.calcSummaryOutput();
                    if (mtt) myMTTOutput.calcSummaryOutput();
                    if (md) myMDOutput.writeSOutputCSV();
                    if (mtt) myMTTOutput.writeSOutputCSV();
                }
            }
        }
    }
    cout << "All tests completed." << endl;
}


int main() {
//    runTest({30}, {40}, {1}, {1}, {1});
//    runTest({30, 50}, {20}, {1, 2, 3}, {1, 2}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
//    runTest({30, 50}, {20}, {1, 2, 3}, {1, 2}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 0, 0, 1, 0, 1);
    //runTest({ 1000 }, { 20,40 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5 }, 0, 1, 0, 1, 0);
//    runTest({1000}, {40}, {3}, {2}, {1, 2, 3, 4, 5}, 0, 1, 0, 1, 0);
//    runTest({ 20 }, { 1000 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
//    runTest({ 100 }, { 20 }, { 2 }, { 1 }, { 1 }, 0, 0, 1, 0, 1);
//    runTest({30}, {20, 60, 100}, {2, 3}, {1, 2}, {1, 2, 3, 4, 5}, 0, 0, 1, 0, 1);
    runTest({50}, {40}, {1, 2, 3}, {1, 2}, {1, 2, 3}, 0, 0, 1, 0, 1);

//    MDP
//    runTest({30}, {60, 100}, {1, 2, 3}, {1, 2}, {1, 2, 3, 4, 5}, 0, 0, 1, 0, 1);
//    runTest({30}, {40}, {1, 2, 3}, {1, 2}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 0, 0, 1, 0, 1);
//    MTTP
    //    runTest({ 50 }, { 20,40,60,80,100 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5,6,7,8,9,10 }, 0, 0, 1, 0, 1);

//    runTest({ 50 }, { 50,100,150,200,250 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
//     runTest({ 10000 }, { 40 }, { 1 }, { 1 }, { 1 }, 0, 1, 0, 1, 0);


    //TEN my_ten(30, 40, 1, 1, 1);
    //Output output = my_ten.findMTT();
    return 0;
}

//