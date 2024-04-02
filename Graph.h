//
// Created by delus on 3/29/24.
//

#ifndef TDSPP_GRAPH_H
#define TDSPP_GRAPH_H

#include <string>
#include <vector>
#include <map>
#include <math.h>

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

    Graph(const int n, const int eT, const int sT = 0, const int start_n = 0, const int end_n = -1);

    void addNode(const int i);

    void addNodes(const int n);

    void addArc(const int start_idx, const int end_idx, vector<double> travel_times);

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
    map<TypeArc, vector<vector<double>>> minttMap; //returns table of min travel time for arc in interval
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
    double TT(const int start_idx, const int end_idx, const double start_t); // TT = Time Travel

    double iTT(const int start_idx, const int end_idx, const double end_t); // iTT = Inverse Time Travel

    double FSP(const vector<double> &start_ts, const vector<double> &end_ts); // Forward shortest path

    vector<double> FSPT(const int start_idx, const double start_t, vector<TypeArc> &arcs);

    vector<double> BSPT(const int end_idx, const double end_t, vector<TypeArc> &arcs);

    double minTT(const int start_idx, const int end_idx, double left_t, double right_t);

    int min_idx_TT(const int start_idx, const int end_idx, const int left_t, const int right_t);

    bool checkFIFO();
};


#endif //TDSPP_GRAPH_H
