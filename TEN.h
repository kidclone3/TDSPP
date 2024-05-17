#ifndef TDSPP_TEN_H
#define TDSPP_TEN_H

#include <vector>
#include <set>
#include <math.h>
#include <iostream>
#include <climits>
#include "Output.h"
#include "Graph.h"

//#include <map>
#include <map>
#include <unordered_map>

//#define map map

using namespace std;

typedef pair<int, int> TypeArc;
typedef pair<int, int> TypeBP;

class TEN {
public:
    //Common
    //Minimum Duration Structures
    struct MDTimedNode { // Minimum Duration Timed Node
        double time = 0;

        MDTimedNode(double t) {
            time = t;
        }
    };

    struct MDTimedNode_compare { // Compare for set data structure
        bool operator()(const MDTimedNode &lhs, const MDTimedNode &rhs) const {
            return lhs.time < rhs.time;
        }
    };

    vector<set<MDTimedNode, MDTimedNode_compare>> timedNodes;

    struct Abspt { // Arc-completed backward shortest path tree
        vector<const MDTimedNode *> nodes;
        vector<double> times;
        double lb = 0;
        double ub = INT_MAX;
        int bpNode = -1; // Breakpoint node
        int bpTime = -1; // Breakpoint time

        Abspt() : nodes(0), times(0) {
            bpNode = 0;
            bpTime = 0;
        }

        Abspt(int i, int t, int n) : nodes(n), times(n) {
            bpNode = i;
            bpTime = t;
        }
    };

    struct Abspt_compare { // Compare for set data structure
        bool operator()(const Abspt &lhs, const Abspt &rhs) const {
            return lhs.times[0] < rhs.times[0];
        }
    };

    struct Abspt_compare_lb { // Compare for set data structure
        bool operator()(const set<Abspt, Abspt_compare>::iterator &lhs,
                        const set<Abspt, Abspt_compare>::iterator &rhs) const {
            return lhs->lb < rhs->lb;
        }
    };

    struct Abspt_compare_ub { // Compare for set data structure
        bool operator()(const set<Abspt, Abspt_compare>::iterator &lhs,
                        const set<Abspt, Abspt_compare>::iterator &rhs) const {
            return lhs->ub < rhs->ub;
        }
    };

    set<Abspt, Abspt_compare> abspts;
    typedef set<Abspt, Abspt_compare>::iterator Abspt_it;
    set<Abspt_it, Abspt_compare_lb> abspt_lbs;
    set<Abspt_it, Abspt_compare_ub> abspt_ubs;

    //Minimum Duration Functions
    void addABSPT(int bpNode, int bpTime);

    void resolveABSPT(Abspt &curr);

    TypeBP findBP(Abspt &curr, int option = 0);

    Output findMD();

    //Minimum Travel Time Structures
    struct Mangrove;

    struct TimedNode {
        double time = 0;
        int nodeID = -1;
        int bpNode = -1;
        int bpTime = -1;
        bool forward = 0;

        TimedNode(int n, double t, int i, int k, bool forward) {
            nodeID = n;
            time = t;
            bpNode = i;
            bpTime = k;
            this->forward = forward;
        }
    };

    struct TimedNodeHash {
        std::size_t operator()(const TimedNode *k) const {
//            return ((hash<int>()(k->nodeID)
//                     ^ (hash<double>()(k->time) << 1)) >> 1)
//                   ^ (hash<int>()(k->bpNode) << 1)
//                   ^ (hash<int>()(k->bpTime) << 1)
//                   ^ (hash<bool>()(k->forward) << 1);
            return reinterpret_cast<std::size_t>(k);
        }
    };

    struct TimedNodeEqual {
        bool operator()(const TimedNode *lhs, const TimedNode *rhs) const {
            return lhs->nodeID == rhs->nodeID
                   && lhs->time == rhs->time
                   && lhs->bpNode == rhs->bpNode
                   && lhs->bpTime == rhs->bpTime
                   && lhs->forward == rhs->forward;
        }
    };

    struct TimedNode_compare {
        bool operator()(const TimedNode &lhs, const TimedNode &rhs) const {
            return lhs.time < rhs.time;
        }
    };

    const TimedNode *origin = nullptr;
    const TimedNode *destination = nullptr;
    typedef set<TimedNode, TimedNode_compare> set_timedNode;
    vector<vector<set_timedNode>> fTimedNodes; //takes ownership of forward nodes, used to add waiting arcs, indexed by [node][bpNode];
    vector<vector<set_timedNode>> bTimedNodes; //takes ownership of backward nodes, used to add waiting arcs, indexed by [node][bpNode];
    struct Mangrove {
        mutable vector<const TimedNode *> f_nodes;
        mutable vector<const TimedNode *> b_nodes;
        mutable vector<double> f_times;
        mutable vector<double> b_times;
        mutable int bpNode = -1;
        int bpTime = -1;
        mutable bool resolved = 0;

        Mangrove(int i, int t, int n) : f_nodes(n), b_nodes(n), f_times(n), b_times(n) {
            bpNode = i;
            bpTime = t;
        }
    };

    struct Mangrove_compare {
        bool operator()(const Mangrove &lhs, const Mangrove &rhs) const {
            return lhs.bpTime < rhs.bpTime;
        }
    };

    vector<set<Mangrove, Mangrove_compare>> mangroves; //takes ownership of mangroves, indexed by [bpNode];
    vector<map<int, const Mangrove *>> mangroveMap; //finds mangrove, indexed by [bpNode][bpTime];
    typedef pair<const TimedNode *, const TimedNode *> TimedArc_ptr;
    map<TimedArc_ptr, double> ttMapLB;
    map<TimedArc_ptr, double> ttMapUB;
    map<const TimedNode *, vector<const TimedNode *>> outMapLB; //internal arcs for LB
    map<const TimedNode *, vector<const TimedNode *>> outMapUB; //internal arcs for UB
    typedef pair<const TimedNode *, double> TnTT;

    struct TnTT_compare {
        bool operator()(const TnTT &lhs, const TnTT &rhs) const {
            return lhs.second > rhs.second;
        }
    };

    typedef vector<const TimedNode *> type_path;

    //Minimum Travel Time Functions
    const Mangrove *addMangrove(int bpNode, int bpTime);

    pair<type_path, double> findLB();

    pair<type_path, double> findUB();

    vector<TypeBP> findBP(const Mangrove *curr, int addmult = 4, int option = 0);

    set<TypeBP> findBP(type_path &path);

    bool isResolved(type_path &path);

    Output findMTT();

    //Enumeration Functions
    Output findEnumMD();

    const Mangrove *addEnumMangrove(int bpNode, int bpTime);

    Output findEnumMTT();

    //Print Functions
    void printCurrentABSPTs();

    void printOptPath(const Abspt_it &abspt, Output &output, bool flag = 0);

    void printOptPath(const Mangrove &mangrove, Output &output, bool flag = 0);

    void printOptMD(const Abspt_it &abspt) {
        cout << "Optimal MD starting time:" << abspt->times[G.startN] << endl;
        cout << "Optimal MD val:" << abspt->lb << endl;
        cout << "Optimal MD breakpoint: (" << abspt->bpNode << ',' << abspt->bpTime << ')' << endl;
    }

    void printPathInfo(pair<type_path, double> &pathLB, pair<type_path, double> &pathUB);

    void printLBPath(type_path path, Output &output, bool flag = 0);

    void printUBPath(type_path path, Output &output, bool flag = 0);

    void printCurrentMangroves();

    void printBPtoAdd(set<TypeBP> &next_bps);

    void printBar() {
        cout << "====================================================================" << endl;
    }

    //Initializer
    Graph G;

    TEN(int n, int eT, int gtype, int tType, int seed, int sT = 0);
};


#endif //TDSPP_TEN_H
