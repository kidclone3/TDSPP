// TDSPP.cpp : This file contains the 'main' function. Program execution begins and ends there.
// TO DO: add more breakpoints associated with minimums (not just THE minimum).


#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>
#include <stack>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <climits>

using namespace std;
typedef pair<int, int> type_arc;
typedef pair<int, int> type_bp;

struct Output {
    bool md = 0;
    bool ddd = 0;
    int bpExplored = 0;
    double runtime = 0;
    double addruntime = 0;
    double spruntime = 0;
    int arctotal = 0;
    int subPathTotal = 1;
    double optVal = 0;
    int iter = 0;

    void writeOutputCSV(string &filename) {
        ofstream myFile(filename, ofstream::app);
        myFile << this->md << ',';
        myFile << this->ddd << ',';
        myFile << this->bpExplored << ',';
        myFile << this->runtime << ',';
        myFile << this->addruntime << ',';
        myFile << this->spruntime << ',';
        myFile << this->arctotal << ',';
        myFile << this->subPathTotal << ',';
        myFile << this->iter << ',';
        myFile << this->optVal << endl;
        myFile.close();
        return;
    }
};

class SummaryOutput {
public:
    bool md;
    int n;
    int T;
    int gtype;
    int ttype;
    int numseed;
    int bpexplored = 0;
    double runtime = 0;
    double enumtime = 0;
    int arctotal = 0;
    int subpathtotal = 0;
    double avgbpexplored = 0;
    int bptotal = 0;
    double percentbp = 0;
    double avgruntime = 0;
    double avgenumtime = 0;
    double percenttime = 0;
    double gavgpercenttime = 1;
    double avgarc = 0;
    double avgsubpath = 0;
    double totaladdruntime = 0;
    double totalspruntime = 0;
    double avgaddruntime = 0;
    double avgspruntime = 0;
    int iter = 0;
    string filename;

    SummaryOutput(string filename, int n, int T, int gtype, int ttype, int numseed, bool md) {
        this->filename = filename;
        this->n = n;
        this->T = T;
        this->gtype = gtype;
        this->ttype = ttype;
        this->numseed = numseed;
        this->md = md;
    }

    void updateSummaryOutput(Output &output) {
        if (output.ddd) {
            bpexplored += output.bpExplored;
            runtime += output.runtime;
            totaladdruntime += output.addruntime;
            totalspruntime += output.spruntime;
            arctotal += output.arctotal;
            subpathtotal += output.subPathTotal;
            gavgpercenttime *= output.runtime;
            iter += output.iter;
        } else {
            bptotal = output.bpExplored;
            enumtime += output.runtime;
            gavgpercenttime /= output.runtime;
        }
        return;
    }

    void calcSummaryOutput() {
        avgarc = (double) arctotal / numseed;
        avgbpexplored = (double) bpexplored / numseed;
        avgenumtime = enumtime / numseed;
        avgruntime = runtime / numseed;
        avgaddruntime = totaladdruntime / numseed;
        avgspruntime = totalspruntime / numseed;
        avgsubpath = (double) subpathtotal / numseed;
        percenttime = 100 * avgruntime / avgenumtime;
        percentbp = 100 * (double) avgbpexplored / bptotal;
        return;
    }

    void writeSOutputCSV() {
        ofstream myFile(filename, ofstream::app);
        myFile << n << ',';
        myFile << T << ',';
        myFile << gtype << ',';
        myFile << ttype << ',';
        myFile << md << ',';
        myFile << avgbpexplored << ',';
        myFile << bptotal << ',';
        myFile << percentbp << ',';
        myFile << pow(gavgpercenttime, 1.0 / numseed) << ',';
        myFile << avgaddruntime << ',';
        myFile << avgspruntime << ',';
        myFile << avgarc << ',';
        myFile << avgsubpath << ',';
        myFile << avgruntime << ',';
        myFile << avgenumtime << ',';
        myFile << percenttime << ',';
        myFile << iter << endl;
        myFile.close();
        return;
    }
};

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
    myFile << "gtype" << ',';
    myFile << "ttype" << ',';
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

    //Initialization
    Graph(const int n, const int eT, const int sT = 0, const int start_n = 0, const int end_n = -1) {
        startT = sT;
        endT = eT;
        this->n = n;
        startN = start_n;
        if (end_n != -1) endN = end_n;
        else endN = n - 1;
        addNodes(n);
    }

    void addNode(const int i) {
        Node newNode(i);
        nodes.push_back(newNode);
    }

    void addNodes(const int n) {
        for (int i = 0; i < n; i++) {
            addNode(i);
        }
    }

    void addArc(const int start_idx, const int end_idx, vector<double> travel_times) {

        /*
        adds arc information to maps
        to get ifloortraveltimes:
            for each integer going forwards, see where it lands, if it overshoots, assign the previous to ifloor
        to get iceiltraveltimes:
            for each integer going backwards, see where it lands, if it undershoots, assign the next to iceil
        to get minttMap:
            

        */
        //Simple Maps
        type_arc arc = {start_idx, end_idx};
        inMap[end_idx].push_back(start_idx);
        outMap[start_idx].push_back(end_idx);
        ttMap[arc] = travel_times;
        //ifloorMap
        int T = endT - startT + 1;
        int curr = 0;
        vector<int> ifloortraveltimes(T);
        for (int i = 0; i < T; i++) {
            while (curr <= T - 1 && curr + travel_times[curr] <= i) {
                curr++;
            }
            ifloortraveltimes[i] = curr - 1;
            //cout << "floor:";
            //cout << i << ',' << travel_times[i] << ',';
            //cout << i << ',' << ifloortraveltimes[i] << endl;
        }
        ittfloorMap[arc] = ifloortraveltimes;
        //iceilMap
        curr = T - 1;
        vector<int> iceiltraveltimes(T);
        for (int i = T - 1; i >= 0; i--) {
            while (curr >= 0 && curr + travel_times[curr] >= i) {
                curr--;
            }
            iceiltraveltimes[i] = curr + 1;
            //cout << "ceil:";
            //cout << i << ',' << travel_times[i] << ',';
            //cout << i << ',' << iceiltraveltimes[i] << endl;
        }
        ittceilMap[arc] = iceiltraveltimes;
        //minttMap and min_idx_ttMap
        vector<vector<double>> mintt(T, vector<double>(T, INT_MAX));
        vector<vector<int>> minindtt(T, vector<int>(T, -1));
        for (int left = 0; left < T; left++) {
            double currmin = travel_times[left];
            int currind = left;
            for (int right = left; right < T; right++) {
                if (currmin > travel_times[right]) {
                    currmin = travel_times[right];
                    currind = right;
                }
                mintt[left][right] = currmin;
                minindtt[left][right] = currind;
            }
        }
        minttMap[arc] = mintt;
        min_idx_ttMap[arc] = minindtt;
    }

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
    map<type_arc, vector<double>> ttMap;
    map<type_arc, vector<int>> ittfloorMap; //returns floor guess for each integer endt
    map<type_arc, vector<int>> ittceilMap; //returns ceil guess for each integer endt
    map<type_arc, vector<vector<double>>> minttMap; //returns table of min travel time for arc in interval
    map<type_arc, vector<vector<int>>> min_idx_ttMap; //returns table of departure time of min travel times for arc in interval
    typedef pair<int, double> type_timednode;

    struct timednode_less {
        bool operator()(const type_timednode &lhs, const type_timednode &rhs) const {
            return lhs.second < rhs.second;
        }
    };

    struct timednode_greater {
        bool operator()(const type_timednode &lhs, const type_timednode &rhs) const {
            return lhs.second > rhs.second;
        }
    };

    //Functions
    double TT(const int start_idx, const int end_idx, const double start_t) { // TT = Time Travel
        /*
        returns travel time of arc from start_idx to end_idx with departure time start_t
        procedure:
        if start_t<startT use travel time at startT
        if start_t>endT use travel time at endT
        otherwise interpolate travel time between floor(start_t) and ceil(start_t)
        */
        pair<int, int> arc = {start_idx, end_idx};
        if (start_t < startT) return ttMap[arc][0];
        else if (start_t > endT) return ttMap[arc][endT - startT];
        else if (start_t == (int) start_t) return ttMap[arc][(int) start_t];
        else {
            double fract_part, int_part;
            fract_part = modf(start_t - startT, &int_part);
            return ttMap[arc][(int) int_part] +
                   fract_part * (ttMap[arc][ceil(start_t) - startT] - ttMap[arc][(int) int_part]);
        }
    }

    double iTT(const int start_idx, const int end_idx, const double end_t) { // iTT = Inverse Time Travel
        /*
        returns travel time of arc from start_idx to end_idx with arrival time end_t
        procedure:
        if end_t<startT+ttMap[arc][0] use travel time at startT
        if end_t>endT+ttMap[arc][endT-startT] use travel time at endT
        otherwise interpolate travel time between floor*(end_t) and ceil*(end_t)
        (separate case for end_t integer)
        where * indicates rounding to nearest endpoint of breakpoints, see below.

        get guess of left and right breakpoint: ittMap[arc][floor(end_t)],ittMap[arc][floor(end_t)]
        refine guesses until we get correct breakpoints leftPrev and rightPrev
        get leftNext and rightNext

        find a such that ln+a(rn-ln)=end_t
        a=(end_t-ln)/(rn-ln)
        return ttMap[arc][lp+a(rp-lp)]=ttMap[arc][lp+(end_t-ln)/(rn-ln)] since rp-lp=1
        

        */
        pair<int, int> arc = {start_idx, end_idx};
        if (end_t < startT + ttMap[arc][0]) return ttMap[arc][0];
        if (end_t > endT + ttMap[arc][endT - startT]) return ttMap[arc][endT - startT];
        int leftPrev, rightPrev;
        if (end_t > endT) {
            leftPrev = ittfloorMap[arc][endT - startT];
            rightPrev = ittceilMap[arc][endT - startT];
        } else {
            leftPrev = ittfloorMap[arc][floor(end_t)];
            rightPrev = ittceilMap[arc][ceil(end_t)];
        }
        if (leftPrev == rightPrev) return ttMap[arc][leftPrev];
        //cout << "before" << leftPrev << ',' << rightPrev << endl;
        while (rightPrev - leftPrev > 1) {
            if (leftPrev + 1 + ttMap[arc][leftPrev + 1] <= end_t) leftPrev++;
            if (rightPrev - 1 + ttMap[arc][rightPrev - 1] > end_t) rightPrev--;
        }
        //cout << "after" << leftPrev << ',' << rightPrev << endl;
        double leftNext = leftPrev + ttMap[arc][leftPrev];
        double rightNext = rightPrev + ttMap[arc][rightPrev];
        //cout << "next" << leftNext << ',' << rightNext << endl;
        return TT(start_idx, end_idx, leftPrev + (end_t - leftNext) / (rightNext - leftNext));
    }

    double FSP(const vector<double> &start_ts, const vector<double> &end_ts) { // Forward shortest path
        vector<double> fspt(n, endT + 1); // Forward shortest path tree
        priority_queue<type_timednode, vector<type_timednode>, timednode_greater> pq;
        pq.push({startN, start_ts[startN]});
        fspt[startN] = start_ts[startN];
        while (!pq.empty()) {
            int i = pq.top().first;
            if (i == endN) break;
            //Care double comparison
            if (pq.top().second != fspt[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (int j: outMap[i]) {
                //cout << i << j << start_ts[i] << end_ts[i] << endl;
                double weight = minTT(i, j, start_ts[i], end_ts[i]);
                if (fspt[j] > fspt[i] + weight) {
                    fspt[j] = fspt[i] + weight;
                    pq.push({j, fspt[j]});
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
        }
        return fspt[endN] - fspt[startN];

    }

    vector<double> FSPT(const int start_idx, const double start_t, vector<type_arc> &arcs) {
        vector<double> fspt(n, endT + 1);
        vector<int> pred(n, -1);
        priority_queue<type_timednode, vector<type_timednode>, timednode_greater> pq;
        pq.push({start_idx, start_t});
        fspt[start_idx] = start_t;
        while (!pq.empty()) {
            int i = pq.top().first;
            //Care double comparison
            if (pq.top().second != fspt[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (int j: outMap[i]) {
                double weight = TT(i, j, fspt[i]);
                if (fspt[j] > fspt[i] + weight) {
                    fspt[j] = fspt[i] + weight;
                    pred[j] = i;
                    pq.push({j, fspt[j]});
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
        }
        for (int j = 0; j < n; j++) {
            if (pred[j] == -1) continue;
            else {
                arcs.push_back({pred[j], j});
            }
        }
        //printf("Times of Nodes in FSPT\n");
        //for (int i = 0; i < n; i++)
        //    printf("%d \t\t %f\n", i, fspt[i]);
        //printf("Arcs in FSPT\n");
        //for (int ind = 0; ind < arcs.size(); ind++)
        //    printf("%d \t\t %d\n", arcs[ind].first, arcs[ind].second);
        return fspt;
    }

    vector<double> BSPT(const int end_idx, const double end_t, vector<type_arc> &arcs) {
        vector<double> bspt(n, startT - 1); // Backward shortest path tree
        vector<int> successor(n, -1);
        priority_queue<type_timednode, vector<type_timednode>, timednode_less> pq;
        pq.push({end_idx, end_t});
        bspt[end_idx] = end_t;
        while (!pq.empty()) {
            int j = pq.top().first;
            //Care double comparison
            if (pq.top().second != bspt[j]) {
                //cout << "(not) popping:" << j << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << j << ',' << bspt[j] << endl;
            pq.pop();
            for (int i: inMap[j]) {
                double weight = iTT(i, j, bspt[j]);
                if (bspt[i] < bspt[j] - weight) {
                    bspt[i] = bspt[j] - weight;
                    successor[i] = j;
                    pq.push({i, bspt[i]});
                    //cout << "pushing:" << i << ',' << bspt[i] << endl;
                }
            }
            //printf("Times of Nodes in BSPT\n");
            //for (int j = 0; j < n; j++)
            //    printf("%d \t\t %f\n", j, bspt[j]);
        }
        for (int i = 0; i < n; i++) {
            if (successor[i] == -1) continue;
            else {
                arcs.push_back({i, successor[i]});
            }
        }
        //printf("Times of Nodes in BSPT\n");
        //for (int j = 0; j < n; j++)
        //    printf("%d \t\t %f\n", j, bspt[j]);
        //printf("Arcs in BSPT\n");
        //for (int ind = 0; ind < arcs.size(); ind++)
        //    printf("%d \t\t %d\n", arcs[ind].first, arcs[ind].second);
        return bspt;
    }

    double minTT(const int start_idx, const int end_idx, double left_t, double right_t) {
        /* Finds departure time of minimum travel time for arc (start_idx,end_idx) in interval [left_t,right_t]*/
        left_t = min((double) endT - startT, max(0.0, left_t - startT));
        right_t = min((double) endT - startT, max(0.0, right_t - startT));
        double ans = min(TT(start_idx, end_idx, left_t), TT(start_idx, end_idx, right_t));
        return min(ans, minttMap[{start_idx, end_idx}][ceil(left_t)][floor(right_t)]);
    }

    int min_idx_TT(const int start_idx, const int end_idx, const int left_t, const int right_t) {
        return min_idx_ttMap[{start_idx, end_idx}][left_t][right_t];
    }

    bool checkFIFO() {
        for (int i = 0; i < n; i++) {
            for (int j: outMap[i]) {
                for (int t = 1; t <= endT - startT; t++) {
                    if (ttMap[{i, j}][t] + 1 < ttMap[{i, j}][t - 1]) {
                        cout << i << ',' << j << endl;
                        cout << t - 1 << ',' << ttMap[{i, j}][t - 1] << endl;
                        cout << t << ',' << ttMap[{i, j}][t] << endl;
                        return 0;
                    }
                }
            }
        }
        return 1;
    }
};

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
    void addABSPT(const int bpNode, const int bpTime) {
        //cout << "Adding ABSPT:(" << bpNode << ',' << bpTime << ')' << endl;
        Abspt new_abspt(bpNode, bpTime, G.n);
        vector<type_arc> arcs = {};
        double endTime = G.FSPT(bpNode, bpTime, arcs)[G.endN];
        arcs = {};
        new_abspt.times = G.BSPT(G.endN, endTime, arcs);
        new_abspt.times[bpNode] = bpTime; //for floating point errors.
        new_abspt.ub = new_abspt.times[G.endN] - new_abspt.times[G.startN];
        for (int i = 0; i < G.n; i++) {
            MDTimedNode timedNode(new_abspt.times[i]);
            auto it = timedNodes[i].insert(timedNode);
            new_abspt.nodes[i] = &(*it.first);
        }
        //Find LB for new ABSPT
        set<Abspt, Abspt_compare>::iterator it = abspts.upper_bound(new_abspt);
        if (it == abspts.end()) {
            new_abspt.lb = new_abspt.ub;
        } else {
            new_abspt.lb = G.FSP(new_abspt.times, it->times);
        }
        //Update LB for prev ABSPT
        if (it == abspts.begin()) {
            cout << "no previous ABSPT, this should only occur during first or second ABSPT" << endl;
        } else {
            it--;
            Abspt prevAbspt = *it;
            abspt_lbs.erase(it);
            abspt_ubs.erase(it);
            abspts.erase(*it);
            prevAbspt.lb = G.FSP(prevAbspt.times, new_abspt.times);
            auto it_bool = abspts.insert(prevAbspt);
            it = it_bool.first;
            //cout << "Updating ABSPT:(" << it->bpNode << ',' << it->bpTime << ')' << endl;
            abspt_lbs.insert(it);
            abspt_ubs.insert(it);
        }
        auto it_bool = abspts.insert(new_abspt);
        it = it_bool.first;
        //cout << "Inserting ABSPT:(" << it->bpNode << ',' << it->bpTime << ')' << endl;
        abspt_lbs.insert(it);
        abspt_ubs.insert(it);
        return;
    }

    void resolveABSPT(Abspt &curr) {
        Abspt currAbspt = curr;
        Abspt_it it = abspts.find(curr);
        abspt_lbs.erase(it);
        abspt_ubs.erase(it);
        abspts.erase(*it);
        currAbspt.lb = currAbspt.ub;
        auto it_bool = abspts.insert(currAbspt);
        it = it_bool.first;
        //cout << "Resolving ABSPT:(" << it->bpNode << ',' << it->bpTime << ')' << endl;
        abspt_lbs.insert(it);
        abspt_ubs.insert(it);
    }

    type_bp findBP(Abspt &curr, int option = 1) {
        auto it = abspts.find(curr);
        it++;
        int bpNode = -1;
        int bpTime = -1;
        if (it == abspts.end()) return {bpNode, bpTime};
        for (int i = 0; i < G.n; i++) {
            if (G.outMap[i].empty()) continue;
            double dleftT = curr.times[i];
            double drightT = it->times[i];
            int leftT, rightT;
            if (dleftT == (int) dleftT) leftT = dleftT + 1;
            else leftT = ceil(dleftT);
            if (drightT == (int) drightT) rightT = drightT - 1;
            else rightT = floor(drightT);
            if (leftT > rightT) continue;
            else {
                bpNode = i;
                if (option == 0) bpTime = G.min_idx_TT(i, G.outMap[i][0], leftT, rightT); //min
                else if (option == 1) bpTime = (leftT + rightT) / 2; //med
                else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
                return {bpNode, bpTime};
            }
        }
        return {bpNode, bpTime};
    }

    Output findMD() {
        int bpExplored = 0;
        int iter = 0;
        auto start = chrono::high_resolution_clock::now();
        addABSPT(G.endN, G.endT);
        addABSPT(G.startN, G.startT);
        iter += 2;
        bpExplored += 2;
        auto lb_it = abspt_lbs.begin();
        auto ub_it = abspt_ubs.begin();
        //cout << "current lower bound=" << (*lb_it)->lb << endl;
        //cout << "current upper bound=" << (*ub_it)->ub << endl;
        //printCurrentABSPTs();
        while ((*lb_it)->lb != (*ub_it)->ub) {
            auto my_abspt = **lb_it;
            type_bp next_bp = findBP(my_abspt);
            if (next_bp.first == -1) {
                resolveABSPT(my_abspt);
            } else {
                addABSPT(next_bp.first, next_bp.second);
                bpExplored++;
            }
            lb_it = abspt_lbs.begin();
            ub_it = abspt_ubs.begin();
            iter++;
            //cout << "current lower bound=" << (*lb_it)->lb << endl;
            //cout << "current upper bound=" << (*ub_it)->ub << endl;
            //printCurrentABSPTs();
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        Output output;
        output.md = 1;
        output.ddd = 1;
        output.bpExplored = bpExplored;
        output.runtime = duration.count();
        output.iter = iter;
        printOptPath(*lb_it, output);
        //printOptMD(*lb_it);
        output.subPathTotal = 1;
        output.optVal = (*lb_it)->lb;
        return output;
    }

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
    const Mangrove *addMangrove(const int bpNode, const int bpTime) {
        //cout << "Adding Mangrove:(" << bpNode << ',' << bpTime << ')' << endl;
        Mangrove new_mangrove(bpNode, bpTime, G.n);
        auto next_it = mangroves[bpNode].upper_bound(new_mangrove);
        const Mangrove *next_man_ptr = (next_it == mangroves[bpNode].end() ? nullptr : &(*next_it));
        const Mangrove *prev_man_ptr = (next_it == mangroves[bpNode].begin() ? nullptr : &(*prev(next_it, 1)));
        ////Forward
        vector<type_arc> f_arcs = {};
        new_mangrove.f_times = G.FSPT(bpNode, bpTime, f_arcs);
        //Add nodes to mangrove and fTimedNodes
        for (int i = 0; i < G.n; i++) {
            if (new_mangrove.f_times[i] <= G.endT) {
                TimedNode f_timedNode(i, new_mangrove.f_times[i], bpNode, bpTime, 1);
                auto it_bool = fTimedNodes[i][bpNode].insert(f_timedNode);
                auto it = it_bool.first;
                auto currTimedNode = &(*it);
                new_mangrove.f_nodes[i] = currTimedNode;
            } else {
                new_mangrove.f_nodes[i] = nullptr;
            }
        }
        ////Backward
        vector<type_arc> b_arcs = {};
        new_mangrove.b_times = G.BSPT(bpNode, bpTime, b_arcs);
        //Add nodes to mangrove and bTimedNodes
        for (int i = 0; i < G.n; i++) {
            if (new_mangrove.b_times[i] >= G.startT) {
                TimedNode b_timedNode(i, new_mangrove.b_times[i], bpNode, bpTime, 0);
                auto it_bool = bTimedNodes[i][bpNode].insert(b_timedNode);
                auto it = it_bool.first;
                auto currTimedNode = &(*it);
                new_mangrove.b_nodes[i] = currTimedNode;
            } else {
                new_mangrove.b_nodes[i] = nullptr;
            }
        }

        //Split into two cases depending on whether added mangrove is resolved(no next mangrove or next mangrove is one unit away)
        //Add entries for outMapLB, outMapUB, ttMapLB, ttMapUB
        if (next_man_ptr == nullptr || next_man_ptr->bpTime == bpTime + 1) {
            //Do everything resolved
            new_mangrove.resolved = 1;
            //Add arcs in FSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
            for (type_arc arc: f_arcs) {
                int i = arc.first, j = arc.second;
                if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
                    const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
                    const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
                    outMapLB[timedNode_i].push_back(timedNode_j);
                    outMapUB[timedNode_i].push_back(timedNode_j);
                    ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
                    ttMapLB[{timedNode_i, timedNode_j}] = ttMapUB[{timedNode_i, timedNode_j}];
                }
            }
            //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
            for (type_arc arc: b_arcs) {
                int i = arc.first, j = arc.second;
                if (new_mangrove.b_times[i] >= G.startT && new_mangrove.b_times[j] >= G.startT) {
                    const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
                    const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
                    outMapLB[timedNode_i].push_back(timedNode_j);
                    outMapUB[timedNode_i].push_back(timedNode_j);
                    ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.b_times[i]);
                    ttMapLB[{timedNode_i, timedNode_j}] = ttMapUB[{timedNode_i, timedNode_j}];
                }
            }
        } else {
            //Do everything normally
            //Add arcs in FSPT to outMapLB and outMapUB
            for (int i = 0; i < G.n; i++) {
                const TimedNode *f_timedNode = new_mangrove.f_nodes[i];
                if (f_timedNode == nullptr) continue;
                for (int j: G.outMap[i]) {
                    if (new_mangrove.f_times[j] <= G.endT) {
                        const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
                        const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
                        outMapLB[timedNode_i].push_back(timedNode_j);
                        ttMapLB[{timedNode_i, timedNode_j}] = G.minTT(i, j, new_mangrove.f_times[i],
                                                                      next_man_ptr->f_times[i]);
                    }
                }
            }
            for (type_arc arc: f_arcs) {
                int i = arc.first, j = arc.second;
                if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
                    const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
                    const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
                    outMapUB[timedNode_i].push_back(timedNode_j);
                    ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
                }
            }
            //Add arcs in BSPT to outMapLB and outMapUB
            for (int i = 0; i < G.n; i++) {
                const TimedNode *b_timedNode = new_mangrove.b_nodes[i];
                if (b_timedNode == nullptr) continue;
                for (int j: G.outMap[i]) {
                    if (new_mangrove.b_times[j] >= G.startT) {
                        const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
                        const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
                        outMapLB[timedNode_i].push_back(timedNode_j);
                        ttMapLB[{timedNode_i, timedNode_j}] = G.minTT(i, j, new_mangrove.b_times[i],
                                                                      next_man_ptr->b_times[i]);
                    }
                }
            }
            for (type_arc arc: b_arcs) {
                int i = arc.first, j = arc.second;
                if (new_mangrove.b_times[i] >= G.startT && new_mangrove.b_times[j] >= G.startT) {
                    const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
                    const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
                    outMapUB[timedNode_i].push_back(timedNode_j);
                    ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.b_times[i]);
                }
            }
        }
        //Update costs for prev mangrove ttMapLB and ttMapUB depending on whether prev mangrove needs to be resolved or not
        if (prev_man_ptr != nullptr) {
            if (prev_man_ptr->bpTime == bpTime - 1) {
                //Resolve
                prev_man_ptr->resolved = 1;
                for (const TimedNode *f_timedNode: prev_man_ptr->f_nodes) {
                    if (f_timedNode == nullptr) continue;
                    outMapLB[f_timedNode] = outMapUB[f_timedNode];
                    for (const TimedNode *next_f_timedNode: outMapUB[f_timedNode]) {
                        TimedArc_ptr timedArc_ptr = {f_timedNode, next_f_timedNode};
                        ttMapLB[timedArc_ptr] = ttMapUB[timedArc_ptr];
                    }
                }
                for (const TimedNode *b_timedNode: prev_man_ptr->b_nodes) {
                    if (b_timedNode == nullptr) continue;
                    outMapLB[b_timedNode] = outMapUB[b_timedNode];
                    for (const TimedNode *next_b_timedNode: outMapUB[b_timedNode]) {
                        TimedArc_ptr timedArc_ptr = {b_timedNode, next_b_timedNode};
                        ttMapLB[timedArc_ptr] = ttMapUB[timedArc_ptr];
                    }
                }
            } else {
                //Update ttMapLB
                for (const TimedNode *f_timedNode: prev_man_ptr->f_nodes) {
                    if (f_timedNode == nullptr) continue;
                    int i = f_timedNode->nodeID;
                    double leftT = prev_man_ptr->f_times[i];
                    for (const TimedNode *next_f_timedNode: outMapLB[f_timedNode]) {
                        int j = next_f_timedNode->nodeID;
                        double rightT = new_mangrove.f_times[j];
                        TimedArc_ptr timedArc_ptr = {f_timedNode, next_f_timedNode};
                        ttMapLB[timedArc_ptr] = G.minTT(i, j, leftT, rightT);
                    }
                }
                for (const TimedNode *b_timedNode: prev_man_ptr->b_nodes) {
                    if (b_timedNode == nullptr) continue;
                    int i = b_timedNode->nodeID;
                    double leftT = prev_man_ptr->b_times[i];
                    for (const TimedNode *next_b_timedNode: outMapLB[b_timedNode]) {
                        int j = next_b_timedNode->nodeID;
                        double rightT = new_mangrove.b_times[j];
                        TimedArc_ptr timedArc_ptr = {b_timedNode, next_b_timedNode};
                        ttMapLB[timedArc_ptr] = G.minTT(i, j, leftT, rightT);
                    }
                }
            }
        }
        auto it_bool = mangroves[bpNode].insert(new_mangrove);
        const Mangrove *man_ptr = &(*it_bool.first);
        mangroveMap[bpNode][bpTime] = man_ptr;
        return man_ptr;
    }

    pair<type_path, double> findLB() {
        map<const TimedNode *, double> dp;
        map<const TimedNode *, const TimedNode *> pred;
        priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
        pq.push({origin, 0});
        dp[origin] = 0;
        while (!pq.empty()) {
            const TimedNode *ptr = pq.top().first;
            //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpNode << ',' << ptr->bpTime << ")" << endl;
            //cout << "DP value: " << pq.top().second << endl;
            if (ptr == destination) break;
            //Care double comparison
            if (pq.top().second > dp[ptr]) {
                pq.pop();
                continue;
            }
            pq.pop();
            for (const TimedNode *next_ptr: outMapLB[ptr]) {
                double weight = ttMapLB[{ptr, next_ptr}];
                auto it = dp.find(next_ptr);
                if (it == dp.end()) dp[next_ptr] = INT_MAX;
                if (dp[next_ptr] > dp[ptr] + weight) {
                    dp[next_ptr] = dp[ptr] + weight;
                    pq.push({next_ptr, dp[next_ptr]});
                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                    //cout << "DP value: " << dp[next_ptr] << endl;
                    pred[next_ptr] = ptr;
                }
            }
            int node = ptr->nodeID;
            int i = ptr->bpNode;
            if (ptr->forward) {
                for (int j = 0; j < G.n; j++) {
                    if (bTimedNodes[node][j].empty()) continue;
                    auto next_it = bTimedNodes[node][j].upper_bound(*ptr);
                    const TimedNode *next_ptr;
                    if (j == i) {
                        //Waiting arc to next copy
                        if (next_it == bTimedNodes[node][j].end()) continue;
                        else next_ptr = &(*next_it);
                    } else {
                        //Waiting arc from i-FSPT to j-BSPT
                        if (next_it == bTimedNodes[node][j].begin()) next_ptr = &(*next_it);
                        else if (next_it == bTimedNodes[node][j].end()) continue;
                        else {
                            auto prev_it = prev(next_it, 1);
                            if (mangroveMap[j][prev_it->bpTime]->resolved) next_ptr = &(*next_it);
                            else next_ptr = &(*prev_it);
                        }
                    }
                    if (next_ptr == nullptr) continue;
                    auto it = dp.find(next_ptr);
                    if (it == dp.end()) dp[next_ptr] = INT_MAX;
                    if (dp[next_ptr] > dp[ptr]) {
                        dp[next_ptr] = dp[ptr];
                        pq.push({next_ptr, dp[next_ptr]});
                        //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                        //cout << "DP value: " << dp[next_ptr] << endl;
                        pred[next_ptr] = ptr;
                    }
                }
            } else {
                if (i == node) {
                    //Waiting arc from BSPT to FSPT through node
                    auto next_it = fTimedNodes[node][i].lower_bound(*ptr);
                    const TimedNode *next_ptr;
                    assert(next_it != fTimedNodes[node][i].end());
                    next_ptr = &(*next_it);
                    auto it = dp.find(next_ptr);
                    if (it == dp.end()) dp[next_ptr] = INT_MAX;
                    if (dp[next_ptr] > dp[ptr]) {
                        dp[next_ptr] = dp[ptr];
                        pq.push({next_ptr, dp[next_ptr]});
                        //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                        //cout << "DP value: " << dp[next_ptr] << endl;
                        pred[next_ptr] = ptr;
                    }
                }
            }
        }
        type_path path;
        const TimedNode *curr = destination;
        while (curr != origin) {
            path.push_back(curr);
            curr = pred[curr];
        }
        path.push_back(curr);
        reverse(path.begin(), path.end());
        return {path, dp[destination] - dp[origin]};
    }

    pair<type_path, double> findUB() {
        map<const TimedNode *, double> dp;
        map<const TimedNode *, const TimedNode *> pred;
        priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
        pq.push({origin, 0});
        dp[origin] = 0;
        while (!pq.empty()) {
            const TimedNode *ptr = pq.top().first;
            //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpNode << ',' << ptr->bpTime << ")" << endl;
            //cout << "DP value: " << pq.top().second << endl;
            if (ptr == destination) break;
            //Care double comparison
            if (pq.top().second > dp[ptr]) {
                pq.pop();
                continue;
            }
            pq.pop();
            for (const TimedNode *next_ptr: outMapUB[ptr]) {
                double weight = ttMapUB[{ptr, next_ptr}];
                auto it = dp.find(next_ptr);
                if (it == dp.end()) dp[next_ptr] = INT_MAX;
                if (dp[next_ptr] > dp[ptr] + weight) {
                    dp[next_ptr] = dp[ptr] + weight;
                    pq.push({next_ptr, dp[next_ptr]});
                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                    //cout << "DP value: " << dp[next_ptr] << endl;
                    pred[next_ptr] = ptr;
                }
            }
            int node = ptr->nodeID;
            int i = ptr->bpNode;
            if (ptr->forward) {
                for (int j = 0; j < G.n; j++) {
                    if (bTimedNodes[node][j].empty()) continue;
                    auto next_it = bTimedNodes[node][j].upper_bound(*ptr);
                    const TimedNode *next_ptr;
                    if (j == i) {
                        //Waiting arc to next copy
                        if (next_it == bTimedNodes[node][j].end()) continue;
                        else next_ptr = &(*next_it);
                    } else {
                        //Waiting arc from i-FSPT to j-BSPT
                        if (next_it == bTimedNodes[node][j].begin()) next_ptr = &(*next_it);
                        else if (next_it == bTimedNodes[node][j].end()) continue;
                        next_ptr = &(*next_it);
                    }
                    if (next_ptr == nullptr) continue;
                    auto it = dp.find(next_ptr);
                    if (it == dp.end()) dp[next_ptr] = INT_MAX;
                    if (dp[next_ptr] > dp[ptr]) {
                        dp[next_ptr] = dp[ptr];
                        pq.push({next_ptr, dp[next_ptr]});
                        //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                        //cout << "DP value: " << dp[next_ptr] << endl;
                        pred[next_ptr] = ptr;
                    }
                }
            } else {
                if (i == node) {
                    //Waiting arc from BSPT to FSPT through node
                    auto next_it = fTimedNodes[node][i].lower_bound(*ptr);
                    const TimedNode *next_ptr;
                    assert(next_it != fTimedNodes[node][i].end());
                    next_ptr = &(*next_it);
                    auto it = dp.find(next_ptr);
                    if (it == dp.end()) dp[next_ptr] = INT_MAX;
                    if (dp[next_ptr] > dp[ptr]) {
                        dp[next_ptr] = dp[ptr];
                        pq.push({next_ptr, dp[next_ptr]});
                        //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
                        //cout << "DP value: " << dp[next_ptr] << endl;
                        pred[next_ptr] = ptr;
                    }
                }
            }
        }
        type_path path;
        const TimedNode *curr = destination;
        while (curr != origin) {
            path.push_back(curr);
            curr = pred[curr];
        }
        path.push_back(curr);
        reverse(path.begin(), path.end());
        return {path, dp[destination] - dp[origin]};
    }

    vector<type_bp> findBP(const Mangrove *curr, int addmult = 4, int option = 0) {
        vector<type_bp> mult_bps;
        int bpNode = curr->bpNode;
        int bpTime = -1;
        auto next_it = mangroves[curr->bpNode].upper_bound(*curr);
        if (next_it == mangroves[curr->bpNode].end()) return mult_bps;
        //Check if any breakpoints in between
        int leftT = curr->bpTime;
        int rightT = next_it->bpTime;
        leftT++, rightT--;
        if (leftT <= rightT) {
            if (addmult == 4) {
                stack<pair<int, int>> st;
                st.push({leftT, rightT});
                while (!st.empty()) {
                    pair<int, int> interval = st.top();
                    st.pop();
                    leftT = interval.first;
                    rightT = interval.second;
                    if (G.outMap[bpNode].empty()) {
                        if (option == 0) bpTime = G.min_idx_TT(G.inMap[bpNode][0], bpNode, leftT, rightT); //min
                        else if (option == 1) bpTime = (leftT + rightT) / 2; //med
                        else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
                    } else {
                        if (option == 0) bpTime = G.min_idx_TT(bpNode, G.outMap[bpNode][0], leftT, rightT); //min
                        else if (option == 1) bpTime = (leftT + rightT) / 2; //med
                        else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
                    }
                    if (bpTime - 1 > leftT) {
                        st.push({leftT, bpTime - 2});
                    }
                    if (bpTime + 1 < rightT) {
                        st.push({bpTime + 2, rightT});
                    }
                    mult_bps.push_back({bpNode, bpTime});
                    if (bpTime - 1 >= leftT) mult_bps.push_back({bpNode, bpTime - 1});
                    if (bpTime + 1 <= rightT) mult_bps.push_back({bpNode, bpTime + 1});
                }
            } else if (addmult == 5) {
                for (int i = leftT; i <= rightT; i++) {
                    mult_bps.push_back({bpNode, i});
                }
            } else {
                if (G.outMap[bpNode].empty()) {
                    if (option == 0) bpTime = G.min_idx_TT(G.inMap[bpNode][0], bpNode, leftT, rightT); //min
                    else if (option == 1) bpTime = (leftT + rightT) / 2; //med
                    else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
                } else {
                    if (option == 0) bpTime = G.min_idx_TT(bpNode, G.outMap[bpNode][0], leftT, rightT); //min
                    else if (option == 1) bpTime = (leftT + rightT) / 2; //med
                    else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
                }
                mult_bps.push_back({bpNode, bpTime});
                if (addmult >= 1 && bpTime - 1 >= leftT) mult_bps.push_back({bpNode, bpTime - 1});
                if (addmult >= 2 && bpTime + 1 <= rightT) mult_bps.push_back({bpNode, bpTime + 1});
                if (addmult >= 3 && leftT < bpTime - 1) mult_bps.push_back({bpNode, leftT});
            }
        }
        return mult_bps;
    }

    set<type_bp> findBP(type_path &path) {
        set<type_bp> bps;
        for (int i = 1; i < path.size(); i++) {
            const TimedNode *timedNode = path[i];
            if (path[i]->bpNode != path[i - 1]->bpNode) {
                //cout << "Not Investigating: (" << timedNode->bpNode << ',' << timedNode->bpTime << ')' << endl;
                continue;
            }
            //cout << "Investigating: (" << timedNode->bpNode << ',' << timedNode->bpTime << ')' << endl;
            vector<type_bp> mult_bp = findBP(mangroveMap[timedNode->bpNode][timedNode->bpTime]);
            for (type_bp bp: mult_bp) {
                bps.insert(bp);
            }
        }
        return bps;
    }

    bool isResolved(type_path &path) {
        bool wait = 1;
        for (int i = 1; i < path.size(); i++) {
            if (path[i]->nodeID == path[i - 1]->nodeID) {
                wait = 1;
            } else if (wait) {
                wait = 0;
                //cout << path[i]->bpNode << path[i]->bpTime << endl;
                if (!mangroveMap[path[i]->bpNode][path[i]->bpTime]->resolved) {
                    return false;
                }
            }
        }
        return true;
    }

    Output findMTT() {
        int bpExplored = 0;
        int iter = 0;
        auto start = chrono::high_resolution_clock::now();
        auto start2 = chrono::high_resolution_clock::now();
        const Mangrove *lastMangrove = addMangrove(G.endN, G.endT);
        auto stop2 = chrono::high_resolution_clock::now();
        auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
        bpExplored++;
        destination = lastMangrove->b_nodes[G.endN];
        for (int i = 0; i < G.n; i++) {
            if (i == G.endN) continue;
            start2 = chrono::high_resolution_clock::now();
            addMangrove(i, floor(lastMangrove->b_times[i]));
            stop2 = chrono::high_resolution_clock::now();
            duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
            bpExplored++;
        }
        start2 = chrono::high_resolution_clock::now();
        const Mangrove *firstMangrove = addMangrove(G.startN, G.startT);
        stop2 = chrono::high_resolution_clock::now();
        duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
        bpExplored++;
        origin = firstMangrove->f_nodes[G.startN];
        for (int i = 0; i < G.n; i++) {
            if (i == G.startN) continue;
            start2 = chrono::high_resolution_clock::now();
            addMangrove(i, ceil(firstMangrove->f_times[i]));
            stop2 = chrono::high_resolution_clock::now();
            duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
            bpExplored++;
        }
        //Add breakpoints associated with minimums or points near minimums
        if (false) {
            for (int i = 0; i < G.n; i++) {
                int j = G.outMap[i][0];
                int leftT = ceil(firstMangrove->f_times[i]) + 1;
                int rightT = floor(lastMangrove->b_times[i]) - 1;
                double reftime = G.minttMap[{i, j}][leftT][rightT];
                queue<vector<int>> q;
                q.push({G.min_idx_ttMap[{i, j}][leftT][rightT], leftT, rightT});
                while (!q.empty()) {
                    vector<int> v = q.front();
                    int t = v[0];
                    leftT = v[1];
                    rightT = v[2];
                    q.pop();
                    start2 = chrono::high_resolution_clock::now();
                    addMangrove(i, t);
                    if (t > leftT) {
                        addMangrove(i, t - 1);
                        bpExplored++;
                    }
                    if (t - 1 > leftT) {
                        addMangrove(i, t - 2);
                        bpExplored++;
                    }
                    if (t < rightT) {
                        addMangrove(i, t + 1);
                        bpExplored++;
                    }
                    if (t + 1 < rightT) {
                        addMangrove(i, t + 2);
                        bpExplored++;
                    }
                    if (t + 2 < rightT) {
                        addMangrove(i, t + 3);
                        bpExplored++;
                    }
                    stop2 = chrono::high_resolution_clock::now();
                    duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
                    bpExplored++;
                    if (leftT < t - 2 && G.minttMap[{i, j}][leftT][t - 3] < 1.03 * reftime) {
                        q.push({G.min_idx_ttMap[{i, j}][leftT][t - 3], leftT, t - 3});
                    }
                    if (rightT > t + 3 && G.minttMap[{i, j}][t + 4][rightT] < 1.03 * reftime) {
                        q.push({G.min_idx_ttMap[{i, j}][t + 4][rightT], t + 4, rightT});
                    }
                }
            }
        }
        //printCurrentMangroves();
        auto start3 = chrono::high_resolution_clock::now();
        pair<type_path, double> pathLB = findLB();
        auto stop3 = chrono::high_resolution_clock::now();
        auto duration3 = chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
        iter++;
        //pair<type_path, double> pathUB = findUB();
        //printPathInfo(pathLB, pathUB);
        //printCurrentMangroves();
        while (!isResolved(pathLB.first)) {
            set<type_bp> next_bps = findBP(pathLB.first);
            //printBPtoAdd(next_bps);
            for (type_bp next_bp: next_bps) {
                start2 = chrono::high_resolution_clock::now();
                addMangrove(next_bp.first, next_bp.second);
                stop2 = chrono::high_resolution_clock::now();
                duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
                bpExplored++;
            }
            start3 = chrono::high_resolution_clock::now();
            pathLB = findLB();
            stop3 = chrono::high_resolution_clock::now();
            //cout << chrono::duration_cast<chrono::milliseconds>(stop3 - start3).count() << endl;
            duration3 += chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
            //pair<type_path, double> pathUB = findUB();
            //printPathInfo(pathLB, pathUB);
            iter++;
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Time taken by adding mangrove: " << duration2.count() << " milliseconds" << endl;
        cout << "Time taken by calc: " << duration3.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << pathLB.second << endl;
        //pair<type_path, double> pathUB = findUB();
        //printPathInfo(pathLB, pathUB);
        Output output;
        output.md = 0;
        output.ddd = 1;
        output.bpExplored = bpExplored;
        output.runtime = duration.count();
        output.spruntime = duration3.count();
        output.addruntime = duration2.count();
        output.iter = iter;
        printLBPath(pathLB.first, output);
        output.optVal = pathLB.second;
        return output;
    }

    //Enumeration Functions
    Output findEnumMD() {
        auto start = chrono::high_resolution_clock::now();
        const Mangrove *opt_mangrove = nullptr;
        double optVal;
        const Mangrove *lastMangrove = addEnumMangrove(G.endN, G.endT);
        opt_mangrove = lastMangrove;
        optVal = lastMangrove->f_times[G.endN] - lastMangrove->b_times[G.startN];
        const Mangrove *firstMangrove = addEnumMangrove(G.startN, G.startT);
        double temp = firstMangrove->f_times[G.endN] - firstMangrove->b_times[G.startN];
        if (temp < optVal) {
            opt_mangrove = firstMangrove;
            optVal = temp;
        }
        for (int i = 0; i < G.n; i++) {
            int leftT, rightT;
            if ((int) lastMangrove->b_times[i] == lastMangrove->b_times[i]) {
                rightT = lastMangrove->b_times[i] - 1;
            } else {
                rightT = floor(lastMangrove->b_times[i]);
            }
            if ((int) firstMangrove->f_times[i] == firstMangrove->f_times[i]) {
                leftT = firstMangrove->f_times[i] + 1;
            } else {
                leftT = ceil(firstMangrove->f_times[i]);
            }
            for (int t = leftT; t <= rightT; t++) {
                const Mangrove *new_mangrove = addEnumMangrove(i, t);
                double temp = new_mangrove->f_times[G.endN] - new_mangrove->b_times[G.startN];
                if (temp < optVal) {
                    opt_mangrove = new_mangrove;
                    optVal = temp;
                }
            }
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << optVal << endl;
        Output output;
        output.md = 1;
        output.ddd = 0;
        output.bpExplored = (G.n - 1) * (G.endT - G.startT - 1) + 2;
        output.runtime = duration.count();
        printOptPath(*opt_mangrove, output);
        output.subPathTotal = 1;
        output.optVal = optVal;
        return output;
    }

    const Mangrove *addEnumMangrove(const int bpNode, const int bpTime) {
        //cout << "Adding Mangrove:(" << bpNode << ',' << bpTime << ')' << endl;
        Mangrove new_mangrove(bpNode, bpTime, G.n);
        ////Forward
        vector<type_arc> f_arcs = {};
        new_mangrove.f_times = G.FSPT(bpNode, bpTime, f_arcs);
        //Add nodes to mangrove and fTimedNodes
        for (int i = 0; i < G.n; i++) {
            if (new_mangrove.f_times[i] <= G.endT) {
                TimedNode f_timedNode(i, new_mangrove.f_times[i], bpNode, bpTime, 1);
                auto it_bool = fTimedNodes[i][bpNode].insert(f_timedNode);
                auto it = it_bool.first;
                auto curr_timedNode = &(*it);
                new_mangrove.f_nodes[i] = curr_timedNode;
            } else {
                new_mangrove.f_nodes[i] = nullptr;
            }
        }
        ////Backward
        vector<type_arc> b_arcs = {};
        new_mangrove.b_times = G.BSPT(bpNode, bpTime, b_arcs);
        //Add nodes to mangrove and bTimedNodes
        for (int i = 0; i < G.n; i++) {
            if (new_mangrove.b_times[i] >= G.startT) {
                TimedNode b_timedNode(i, new_mangrove.b_times[i], bpNode, bpTime, 0);
                auto it_bool = bTimedNodes[i][bpNode].insert(b_timedNode);
                auto it = it_bool.first;
                auto curr_timedNode = &(*it);
                new_mangrove.b_nodes[i] = curr_timedNode;
            } else {
                new_mangrove.b_nodes[i] = nullptr;
            }
        }
        //Do everything resolved
        new_mangrove.resolved = 1;
        //Add arcs in FSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
        for (type_arc arc: f_arcs) {
            int i = arc.first, j = arc.second;
            if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
                const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
                const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
                outMapUB[timedNode_i].push_back(timedNode_j);
                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
            }
        }
        //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
        for (type_arc arc: b_arcs) {
            int i = arc.first, j = arc.second;
            if (new_mangrove.b_times[i] >= G.startT && new_mangrove.b_times[j] >= G.startT) {
                const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
                const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
                outMapUB[timedNode_i].push_back(timedNode_j);
                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.b_times[i]);
            }
        }
        auto it_bool = mangroves[bpNode].insert(new_mangrove);
        const Mangrove *man_ptr = &(*it_bool.first);
        mangroveMap[bpNode][bpTime] = man_ptr;
        return man_ptr;
    }

    Output findEnumMTT() {
        auto start = chrono::high_resolution_clock::now();
        const Mangrove *lastMangrove = addEnumMangrove(G.endN, G.endT);
        destination = lastMangrove->b_nodes[G.endN];
        const Mangrove *firstMangrove = addEnumMangrove(G.startN, G.startT);
        origin = firstMangrove->f_nodes[G.startN];
        for (int i = 0; i < G.n; i++) {
            int leftT, rightT;
            if ((int) lastMangrove->b_times[i] == lastMangrove->b_times[i]) {
                rightT = lastMangrove->b_times[i] - 1;
            } else {
                rightT = floor(lastMangrove->b_times[i]);
            }
            if ((int) firstMangrove->f_times[i] == firstMangrove->f_times[i]) {
                leftT = firstMangrove->f_times[i] + 1;
            } else {
                leftT = ceil(firstMangrove->f_times[i]);
            }
            for (int t = leftT; t <= rightT; t++) {
                const Mangrove *new_mangrove = addEnumMangrove(i, t);
            }
        }
        pair<type_path, double> pathOpt = findUB();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << pathOpt.second << endl;
        Output output;
        output.md = 0;
        output.ddd = 0;
        output.bpExplored = (G.n - 1) * (G.endT - G.startT - 1) + 2;
        output.runtime = duration.count();
        printUBPath(pathOpt.first, output);
        output.optVal = pathOpt.second;
        return output;
    };

    //Print Functions
    void printCurrentABSPTs() {
        printBar();
        cout << "Printing Current ABSPTs:" << endl;
        for (auto abspt: abspts) {
            cout << "(" << abspt.bpNode << ',' << abspt.bpTime << "), resolved = " << (bool) (abspt.lb == abspt.ub);
            cout << ", lb = " << abspt.lb << ", ub = " << abspt.ub << endl;
        }
        printBar();
    }

    void printOptPath(const Abspt_it &abspt, Output &output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing Opt Path: " << endl;
        output.arctotal = 0;
        vector<type_arc> arcs = {};
        vector<double> times = G.FSPT(G.startN, abspt->times[G.startN], arcs);
        int curr = G.endN;
        stack<pair<int, double>> st;
        while (curr != G.startN) {
            for (type_arc arc: arcs) {
                if (curr == arc.second) {
                    st.push({curr, times[curr]});
                    curr = arc.first;
                    break;
                }
            }
        }
        if (flag) cout << "(" << curr << ',' << times[curr] << ")" << endl;
        while (!st.empty()) {
            auto curr_times = st.top();
            if (flag) cout << "(" << curr_times.first << ',' << curr_times.second << ")" << endl;
            output.arctotal++;
            st.pop();
        }
        if (flag) printBar();
    }

    void printOptPath(const Mangrove &mangrove, Output &output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing Opt Path: " << endl;
        output.arctotal = 0;
        vector<type_arc> arcs = {};
        vector<double> times = G.FSPT(G.startN, mangrove.b_times[G.startN], arcs);
        int curr = G.endN;
        stack<pair<int, double>> st;
        while (curr != G.startN) {
            for (type_arc arc: arcs) {
                if (curr == arc.second) {
                    st.push({curr, times[curr]});
                    curr = arc.first;
                    break;
                }
            }
        }
        if (flag) cout << "(" << curr << ',' << times[curr] << ")" << endl;
        while (!st.empty()) {
            auto curr_times = st.top();
            if (flag) cout << "(" << curr_times.first << ',' << curr_times.second << ")" << endl;
            output.arctotal++;
            st.pop();
        }
        if (flag) printBar();
    }

    void printOptMD(const Abspt_it &abspt) {
        cout << "Optimal MD starting time:" << abspt->times[G.startN] << endl;
        cout << "Optimal MD val:" << abspt->lb << endl;
        cout << "Optimal MD breakpoint: (" << abspt->bpNode << ',' << abspt->bpTime << ')' << endl;
    }

    void printPathInfo(pair<type_path, double> &pathLB, pair<type_path, double> &pathUB) {
        Output output;
        printLBPath(pathLB.first, output, 1);
        cout << "current lower bound = " << pathLB.second << endl;
        printUBPath(pathUB.first, output, 1);
        cout << "current upper bound = " << pathUB.second << endl;
    }

    void printLBPath(type_path path, Output &output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing LB Path: " << endl;
        output.arctotal = 0;
        output.subPathTotal = 0;
        bool waiting = 1;
        for (int i = 1; i < path.size(); i++) {
            const TimedNode *node = path[i];
            if (path[i - 1]->nodeID != path[i]->nodeID) {
                TimedArc_ptr arc = {path[i - 1], path[i]};
                if (flag) cout << "arc cost: " << ttMapLB[arc] << endl;
                double UB_cost = (path[i - 1]->nodeID == path[i]->nodeID) ? 0 : path[i]->time - path[i - 1]->time;
                if (flag) cout << "UB arc cost: " << UB_cost << endl;
                output.arctotal++;
                if (flag)
                    cout << "at: (" << node->nodeID << ',' << node->time << "), in mangrove: (" << node->bpNode << ','
                         << node->bpTime << ")" << endl;
                if (waiting) {
                    waiting = 0;
                    output.subPathTotal++;
                }
            } else {
                waiting = 1;
            }
        }
        if (flag) printBar();
    }

    void printUBPath(type_path path, Output &output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing UB Path: " << endl;
        output.arctotal = 0;
        output.subPathTotal = 0;
        bool waiting = 1;
        for (int i = 1; i < path.size(); i++) {
            const TimedNode *node = path[i];
            if (path[i - 1]->nodeID != path[i]->nodeID) {
                TimedArc_ptr arc = {path[i - 1], path[i]};
                if (flag) cout << "arc cost: " << ttMapUB[arc] << endl;
                output.arctotal++;
                if (flag)
                    cout << "at: (" << node->nodeID << ',' << node->time << "), in mangrove: (" << node->bpNode << ','
                         << node->bpTime << ")" << endl;
                if (waiting) {
                    waiting = 0;
                    output.subPathTotal++;
                }
            } else {
                waiting = 1;
            }
        }
        if (flag) printBar();
    }

    void printCurrentMangroves() {
        printBar();
        cout << "Printing Current Mangroves:" << endl;
        for (auto mangrove_set: mangroves) {
            for (auto mangrove: mangrove_set) {
                cout << "(" << mangrove.bpNode << ',' << mangrove.bpTime << "), resolved = " << mangrove.resolved
                     << endl;
            }
        }
        printBar();
    }

    void printBPtoAdd(set<type_bp> &next_bps) {
        printBar();
        cout << "Printing BP to Add:" << endl;
        for (type_bp bp: next_bps) {
            cout << "(" << bp.first << ',' << bp.second << ")" << endl;
        }
        printBar();
    }

    void printBar() {
        cout << "====================================================================" << endl;
    }

    //Initializer
    Graph G;

    TEN(const int n, const int eT, const int gtype, const int tType, const int seed, const int sT = 0) : timedNodes(n),
                                                                                                         fTimedNodes(n,
                                                                                                                     vector<set_timedNode>(
                                                                                                                             n)),
                                                                                                         bTimedNodes(n,
                                                                                                                     vector<set_timedNode>(
                                                                                                                             n)),
                                                                                                         mangroves(n),
                                                                                                         mangroveMap(
                                                                                                                 n),
                                                                                                         G(n, eT, sT) {
        /*number of nodes, end time, graph type, travel time type, seed, start time (optional)*/
        string filename =
                "Data/n" + to_string(n) + "T" + to_string(eT) + "gt" + to_string(gtype) + "tt" + to_string(tType) +
                "s" + to_string(seed) + ".csv";
        cout << "opening:" << filename << endl;
        std::ifstream nodeFile(filename.c_str());
        if (!nodeFile.is_open()) {
            cerr << "error opening file" << endl;
        }
        string headers;
        getline(nodeFile, headers);
        string line;
        string startID;
        string endID;
        string time;
        while (getline(nodeFile, startID, ',')) {
            vector<double> times;
            getline(nodeFile, endID, ',');
            for (int i = sT; i < eT; i++) {
                getline(nodeFile, time, ',');
                times.push_back(stod(time, nullptr));
            }
            getline(nodeFile, time);
            times.push_back(stod(time, nullptr));
            G.addArc(stoi(startID, nullptr), stoi(endID, nullptr), times);
            //cout << "i:" << stoi(startID, nullptr) << "j:" << stoi(endID, nullptr) << "t:" << times[eT-sT] << endl;
        }
        if (G.checkFIFO()) {
            cout << "read successful!" << endl;
        } else {
            cerr << "Error: FIFO violated" << endl;
            cin.get();
        }
        return;
    }
};

void runTest(vector<int> ns, vector<int> endTs, vector<int> gTypes, vector<int> tTypes, vector<int> seeds,
             const int common_startTime = 0, bool md = 1, bool mtt = 1, bool mdEnum = 1, bool mttEnum = 1) {
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
    string md_filename = "Results/experimentsMDMED.csv";
    string mtt_filename = "Results/MTTBPMMED.csv";
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
                    string filename = "Results/experimentsMDMEDlog.csv";
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
}


int main() {
    runTest({30}, {40}, {1}, {1}, {1});
//    runTest({30, 50}, {20}, {1, 2, 3}, {1, 2}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    //runTest({ 1000 }, { 20,40 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5 }, 0, 1, 0, 1, 0);
//    runTest({1000}, {40}, {3}, {2}, {1, 2, 3, 4, 5}, 0, 1, 0, 1, 0);
//    runTest({ 20 }, { 1000 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
//    runTest({ 100 }, { 20 }, { 2 }, { 1 }, { 1 }, 0, 0, 1, 0, 1);
//    runTest({30}, {20, 60, 100}, {2, 3}, {1, 2}, {1, 2, 3, 4, 5}, 0, 0, 1, 0, 1);
    //runTest({ 50 }, { 40 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 0);
    //runTest({ 50 }, { 20,40,60,80,100 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5,6,7,8,9,10 }, 0, 1, 0, 1, 0);
    //runTest({ 60 }, { 50,100,150,200,250 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
    // runTest({ 10000 }, { 40 }, { 1 }, { 1 }, { 1 }, 0, 1, 0, 1, 0);


    //TEN my_ten(30, 40, 1, 1, 1);
    //Output output = my_ten.findMTT();
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
