#include <cassert>
#include <algorithm>
#include <stack>
#include <fstream>
#include <queue>
#include <random>

#include <chrono>
#include "TEN.h"

//Minimum Duration Functions
void TEN::addABSPT(int bpNode, int bpTime) {
    //cout << "Adding ABSPT:(" << bpNode << ',' << bpTime << ')' << endl;
    Abspt new_abspt(bpNode, bpTime, G.n);
    vector<TypeArc> arcs = {};
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

void TEN::resolveABSPT(Abspt &curr) {
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

TypeBP TEN::findBP(Abspt &curr, int option) {
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

Output TEN::findMD() {
    // Initialize breakpoint exploration and iteration counters
    int bpExplored = 0;
    int iter = 0;

    // Start the timer
    auto start = chrono::high_resolution_clock::now();

    // Add the end and start nodes to the ABSPT
    addABSPT(G.endN, G.endT);
    addABSPT(G.startN, G.startT);
    iter += 2;
    bpExplored += 2;

    // Get iterators to the lower and upper bounds
    auto lb_it = abspt_lbs.begin();
    auto ub_it = abspt_ubs.begin();
    //cout << "current lower bound=" << (*lb_it)->lb << endl;
    //cout << "current upper bound=" << (*ub_it)->ub << endl;
    //printCurrentABSPTs();

    // Loop until the lower and upper bounds are equal
    while ((*lb_it)->lb != (*ub_it)->ub) {
        auto my_abspt = **lb_it;

        // Find the next break point
        TypeBP next_bp = findBP(my_abspt);

        // If no breakpoint is found, resolve the ABSPT
        // Otherwise, add the new ABSPT
        if (next_bp.first == -1) {
            resolveABSPT(my_abspt);
        } else {
            addABSPT(next_bp.first, next_bp.second);
            bpExplored++;
        }
//            Update iterators
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

// Prepare output
    Output output;
    output.md = 1;
    output.ddd = 1;
    output.bpExplored = bpExplored;
    output.runtime = duration.count();
    output.iter = iter;

    printOptPath(*lb_it, output);
    //printOptMD(*lb_it);
//        Set the remaining output parameters
    output.subPathTotal = 1;
    output.optVal = (*lb_it)->lb;
    return output;
}

//Enumeration Functions
Output TEN::findEnumMD() {
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

const TEN::Mangrove *TEN::addEnumMangrove(int bpNode, int bpTime) {
    //cout << "Adding Mangrove:(" << bpNode << ',' << bpTime << ')' << endl;
    Mangrove new_mangrove(bpNode, bpTime, G.n);
    ////Forward
    vector<TypeArc> f_arcs = {};
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
    vector<TypeArc> b_arcs = {};
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
    for (TypeArc arc: f_arcs) {
        int i = arc.first, j = arc.second;
        if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
            const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
            const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
            outMapUB[timedNode_i].push_back(timedNode_j);
            ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
        }
    }
    //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
    for (TypeArc arc: b_arcs) {
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

//Output TEN::findEnumMTT() {
//    auto start = chrono::high_resolution_clock::now();
//    const Mangrove *lastMangrove = addEnumMangrove(G.endN, G.endT);
//    destination = lastMangrove->b_nodes[G.endN];
//    const Mangrove *firstMangrove = addEnumMangrove(G.startN, G.startT);
//    origin = firstMangrove->f_nodes[G.startN];
//    for (int i = 0; i < G.n; i++) {
//        int leftT, rightT;
//        if ((int) lastMangrove->b_times[i] == lastMangrove->b_times[i]) {
//            rightT = lastMangrove->b_times[i] - 1;
//        } else {
//            rightT = floor(lastMangrove->b_times[i]);
//        }
//        if ((int) firstMangrove->f_times[i] == firstMangrove->f_times[i]) {
//            leftT = firstMangrove->f_times[i] + 1;
//        } else {
//            leftT = ceil(firstMangrove->f_times[i]);
//        }
//        for (int t = leftT; t <= rightT; t++) {
//            const Mangrove *new_mangrove = addEnumMangrove(i, t);
//        }
//    }
//    pair<type_path, double> pathOpt = findUB();
//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
//    cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
//    cout << "Optimal Solution is: " << pathOpt.second << endl;
//    Output output;
//    output.md = 0;
//    output.ddd = 0;
//    output.bpExplored = (G.n - 1) * (G.endT - G.startT - 1) + 2;
//    output.runtime = duration.count();
//    printUBPath(pathOpt.first, output);
//    output.optVal = pathOpt.second;
//    return output;
//};

//Print Functions
void TEN::printCurrentABSPTs() {
    printBar();
    cout << "Printing Current ABSPTs:" << endl;
    for (auto abspt: abspts) {
        cout << "(" << abspt.bpNode << ',' << abspt.bpTime << "), resolved = " << (bool) (abspt.lb == abspt.ub);
        cout << ", lb = " << abspt.lb << ", ub = " << abspt.ub << endl;
    }
    printBar();
}

void TEN::printOptPath(const Abspt_it &abspt, Output &output, bool flag) {
    if (flag) printBar();
    if (flag) cout << "Printing Opt Path: " << endl;
    output.arc_total = 0;
    vector<TypeArc> arcs = {};
    vector<double> times = G.FSPT(G.startN, abspt->times[G.startN], arcs);
    int curr = G.endN;
    stack<pair<int, double>> st;
    while (curr != G.startN) {
        for (TypeArc arc: arcs) {
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
        output.arc_total++;
        st.pop();
    }
    if (flag) printBar();
}

void TEN::printOptPath(const Mangrove &mangrove, Output &output, bool flag) {
    if (flag) printBar();
    if (flag) cout << "Printing Opt Path: " << endl;
    output.arc_total = 0;
    vector<TypeArc> arcs = {};
    vector<double> times = G.FSPT(G.startN, mangrove.b_times[G.startN], arcs);
    int curr = G.endN;
    stack<pair<int, double>> st;
    while (curr != G.startN) {
        for (TypeArc arc: arcs) {
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
        output.arc_total++;
        st.pop();
    }
    if (flag) printBar();
}


void TEN::printPathInfo(pair<type_path, double> &pathLB, pair<type_path, double> &pathUB) {
    Output output;
    printLBPath(pathLB.first, output, 1);
    cout << "current lower bound = " << pathLB.second << endl;
    printUBPath(pathUB.first, output, 1);
    cout << "current upper bound = " << pathUB.second << endl;
}

void TEN::printLBPath(type_path path, Output &output, bool flag) {
    if (flag) printBar();
    if (flag) cout << "Printing LB Path: " << endl;
    output.arc_total = 0;
    output.subPathTotal = 0;
    bool waiting = 1;
    for (int i = 1; i < path.size(); i++) {
        const TimedNode *node = path[i];
        if (path[i - 1]->nodeID != path[i]->nodeID) {
            TimedArc_ptr arc = {path[i - 1], path[i]};
            if (flag) cout << "arc cost: " << ttMapLB[arc] << endl;
            double UB_cost = (path[i - 1]->nodeID == path[i]->nodeID) ? 0 : path[i]->time - path[i - 1]->time;
            if (flag) cout << "UB arc cost: " << UB_cost << endl;
            output.arc_total++;
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

void TEN::printUBPath(type_path path, Output &output, bool flag) {
    if (flag) printBar();
    if (flag) cout << "Printing UB Path: " << endl;
    output.arc_total = 0;
    output.subPathTotal = 0;
    bool waiting = 1;
    for (int i = 1; i < path.size(); i++) {
        const TimedNode *node = path[i];
        if (path[i - 1]->nodeID != path[i]->nodeID) {
            TimedArc_ptr arc = {path[i - 1], path[i]};
            if (flag) cout << "arc cost: " << ttMapUB[arc] << endl;
            output.arc_total++;
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

void TEN::printCurrentMangroves() {
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

void TEN::printBPtoAdd(set<TypeBP> &next_bps) {
    printBar();
    cout << "Printing BP to Add:" << endl;
    for (TypeBP bp: next_bps) {
        cout << "(" << bp.first << ',' << bp.second << ")" << endl;
    }
    printBar();
}


//Initializer
TEN::TEN(int n, int eT, int gtype, int tType, int seed, int sT) : timedNodes(n),
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
