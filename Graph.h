//
// Created by delus on 3/29/24.
//

#ifndef TDSPP_GRAPH_H
#define TDSPP_GRAPH_H

#include <string>
#include <vector>
#include <math.h>
#include <climits>
#include <queue>
#include <iostream>

#include <map>

//#define map map

using namespace std;

typedef pair<int, int> TypeArc;
typedef pair<int, int> TypeBP;

class Graph {

public:
    //Node
    struct Node {
        int nodeID = 0;
        string nodeName;

        Node(int ID, string name = "") {
            nodeID = ID;
            if (name.empty()) name = to_string(ID);
            nodeName = name;
        }
    };

    Graph(int n, int eT, int sT = 0, int start_n = 0, int end_n = -1);

    void addNode(int i);

    void addNodes(int n);

    void addArc(int start_idx, int end_idx, vector<double> travel_times);

    //Properties
    int startT = 0;
    int endT = 0;
    int n = 0;
    int startN = 0;
    int endN = 0;
    vector<Node> nodes;
    //map<int, Node*> indMap;
    map<int, vector<int>> inMap;
    map<int, vector<int>> outMap;
    map<TypeArc, vector<double>> ttMap;
    map<TypeArc, vector<int>> ittfloorMap; //returns floor guess for each integer endt
    map<TypeArc, vector<int>> ittceilMap; //returns ceil guess for each integer endt
    map<TypeArc, vector<vector<double>>> minttMap; //returns table of min travel time for arc in interval aka UTT
    map<TypeArc, vector<vector<int>>> min_idx_ttMap; //returns table of departure time of min travel times for arc in interval
    typedef pair<int, double> type_timedNode;

    struct timedNode_less {
        bool operator()(const type_timedNode &lhs, const type_timedNode &rhs) const {
            return lhs.second < rhs.second;
        }
    };

    struct timedNode_greater {
        bool operator()(const type_timedNode &lhs, const type_timedNode &rhs) const {
            return lhs.second > rhs.second;
        }
    };

    //Functions
    double TT(int start_idx, int end_idx, double start_t); // TT = Time Travel

    double iTT(int start_idx, int end_idx, double end_t); // iTT = Inverse Time Travel

    double FSP(const vector<double> &start_ts, const vector<double> &end_ts); // Forward shortest path

    vector<double> FSPT(int start_idx, double start_t, vector<TypeArc> &arcs);

    vector<double> BSPT(int end_idx, double end_t, vector<TypeArc> &arcs);

    double minTT(int start_idx, int end_idx, double left_t, double right_t);

    int min_idx_TT(int start_idx, int end_idx, int left_t, int right_t);

    bool checkFIFO();
};


#endif //TDSPP_GRAPH_H
