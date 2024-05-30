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


////Minimum Travel Time Functions
//const TEN::Mangrove *TEN::addMangrove(int bpNode, int bpTime) {
//    //cout << "Adding Mangrove:(" << bpNode << ',' << bpTime << ')' << endl;
//    Mangrove new_mangrove(bpNode, bpTime, G.n);
//    auto next_it = mangroves[bpNode].upper_bound(new_mangrove);
//    const Mangrove *next_man_ptr = (next_it == mangroves[bpNode].end() ? nullptr : &(*next_it));
//    const Mangrove *prev_man_ptr = (next_it == mangroves[bpNode].begin() ? nullptr : &(*prev(next_it, 1)));
//    ////Forward
//    vector<TypeArc> f_arcs = {};
//    new_mangrove.f_times = G.FSPT(bpNode, bpTime, f_arcs);
//    //Add nodes to mangrove and fTimedNodes
//    for (int i = 0; i < G.n; i++) {
//        if (new_mangrove.f_times[i] <= G.endT) {
//            TimedNode f_timedNode(i, new_mangrove.f_times[i], bpNode, bpTime, 1);
//            auto it_bool = fTimedNodes[i][bpNode].insert(f_timedNode);
//            auto it = it_bool.first;
//            auto currTimedNode = &(*it);
//            new_mangrove.f_nodes[i] = currTimedNode;
//        } else {
//            new_mangrove.f_nodes[i] = nullptr;
//        }
//    }
//    ////Backward
//    vector<TypeArc> b_arcs = {};
//    new_mangrove.b_times = G.BSPT(bpNode, bpTime, b_arcs);
//    //Add nodes to mangrove and bTimedNodes
//    for (int i = 0; i < G.n; i++) {
//        if (new_mangrove.b_times[i] >= G.startT) {
//            TimedNode b_timedNode(i, new_mangrove.b_times[i], bpNode, bpTime, 0);
//            auto it_bool = bTimedNodes[i][bpNode].insert(b_timedNode);
//            auto it = it_bool.first;
//            auto currTimedNode = &(*it);
//            new_mangrove.b_nodes[i] = currTimedNode;
//        } else {
//            new_mangrove.b_nodes[i] = nullptr;
//        }
//    }
//
//    //Split into two cases depending on whether added mangrove is resolved(no next mangrove or next mangrove is one unit away)
//    //Add entries for outMapLB, outMapUB, ttMapLB, ttMapUB
//    if (next_man_ptr == nullptr || next_man_ptr->bpTime == bpTime + 1) {
//        //Do everything resolved
//        new_mangrove.resolved = 1;
//        //Add arcs in FSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
//        for (TypeArc arc: f_arcs) {
//            int i = arc.first, j = arc.second;
//            if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
//                const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
//                const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
//                outMapLB[timedNode_i].push_back(timedNode_j);
//                outMapUB[timedNode_i].push_back(timedNode_j);
//                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
//                ttMapLB[{timedNode_i, timedNode_j}] = ttMapUB[{timedNode_i, timedNode_j}];
//            }
//        }
//        //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
//        for (TypeArc arc: b_arcs) {
//            int i = arc.first, j = arc.second;
//            if (new_mangrove.b_times[i] >= G.startT && new_mangrove.b_times[j] >= G.startT) {
//                const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
//                const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
//                outMapLB[timedNode_i].push_back(timedNode_j);
//                outMapUB[timedNode_i].push_back(timedNode_j);
//                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.b_times[i]);
//                ttMapLB[{timedNode_i, timedNode_j}] = ttMapUB[{timedNode_i, timedNode_j}];
//            }
//        }
//    } else {
//        //Do everything normally
//        //Add arcs in FSPT to outMapLB and outMapUB
//        for (int i = 0; i < G.n; i++) {
//            const TimedNode *f_timedNode = new_mangrove.f_nodes[i];
//            if (f_timedNode == nullptr) continue;
//            for (int j: G.outMap[i]) {
//                if (new_mangrove.f_times[j] <= G.endT) {
//                    const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
//                    const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
//                    outMapLB[timedNode_i].push_back(timedNode_j);
//                    ttMapLB[{timedNode_i, timedNode_j}] = G.minTT(i, j, new_mangrove.f_times[i],
//                                                                  next_man_ptr->f_times[i]);
//                }
//            }
//        }
//        for (TypeArc arc: f_arcs) {
//            int i = arc.first, j = arc.second;
//            if (new_mangrove.f_times[i] <= G.endT && new_mangrove.f_times[j] <= G.endT) {
//                const TimedNode *timedNode_i = new_mangrove.f_nodes[i];
//                const TimedNode *timedNode_j = new_mangrove.f_nodes[j];
//                outMapUB[timedNode_i].push_back(timedNode_j);
//                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.f_times[i]);
//            }
//        }
//        //Add arcs in BSPT to outMapLB and outMapUB
//        for (int i = 0; i < G.n; i++) {
//            const TimedNode *b_timedNode = new_mangrove.b_nodes[i];
//            if (b_timedNode == nullptr) continue;
//            for (int j: G.outMap[i]) {
//                if (new_mangrove.b_times[j] >= G.startT) {
//                    const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
//                    const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
//                    outMapLB[timedNode_i].push_back(timedNode_j);
//                    ttMapLB[{timedNode_i, timedNode_j}] = G.minTT(i, j, new_mangrove.b_times[i],
//                                                                  next_man_ptr->b_times[i]);
//                }
//            }
//        }
//        for (TypeArc arc: b_arcs) {
//            int i = arc.first, j = arc.second;
//            if (new_mangrove.b_times[i] >= G.startT && new_mangrove.b_times[j] >= G.startT) {
//                const TimedNode *timedNode_i = new_mangrove.b_nodes[i];
//                const TimedNode *timedNode_j = new_mangrove.b_nodes[j];
//                outMapUB[timedNode_i].push_back(timedNode_j);
//                ttMapUB[{timedNode_i, timedNode_j}] = G.TT(i, j, new_mangrove.b_times[i]);
//            }
//        }
//    }
//    //Update costs for prev mangrove ttMapLB and ttMapUB depending on whether prev mangrove needs to be resolved or not
//    if (prev_man_ptr != nullptr) {
//        if (prev_man_ptr->bpTime == bpTime - 1) {
//            //Resolve
//            prev_man_ptr->resolved = 1;
//            for (const TimedNode *f_timedNode: prev_man_ptr->f_nodes) {
//                if (f_timedNode == nullptr) continue;
//                outMapLB[f_timedNode] = outMapUB[f_timedNode];
//                for (const TimedNode *next_f_timedNode: outMapUB[f_timedNode]) {
//                    TimedArc_ptr timedArc_ptr = {f_timedNode, next_f_timedNode};
//                    ttMapLB[timedArc_ptr] = ttMapUB[timedArc_ptr];
//                }
//            }
//            for (const TimedNode *b_timedNode: prev_man_ptr->b_nodes) {
//                if (b_timedNode == nullptr) continue;
//                outMapLB[b_timedNode] = outMapUB[b_timedNode];
//                for (const TimedNode *next_b_timedNode: outMapUB[b_timedNode]) {
//                    TimedArc_ptr timedArc_ptr = {b_timedNode, next_b_timedNode};
//                    ttMapLB[timedArc_ptr] = ttMapUB[timedArc_ptr];
//                }
//            }
//        } else {
//            //Update ttMapLB
//            for (const TimedNode *f_timedNode: prev_man_ptr->f_nodes) {
//                if (f_timedNode == nullptr) continue;
//                int i = f_timedNode->nodeID;
//                double leftT = prev_man_ptr->f_times[i];
//                for (const TimedNode *next_f_timedNode: outMapLB[f_timedNode]) {
//                    int j = next_f_timedNode->nodeID;
//                    double rightT = new_mangrove.f_times[j];
//                    TimedArc_ptr timedArc_ptr = {f_timedNode, next_f_timedNode};
//                    ttMapLB[timedArc_ptr] = G.minTT(i, j, leftT, rightT);
//                }
//            }
//            for (const TimedNode *b_timedNode: prev_man_ptr->b_nodes) {
//                if (b_timedNode == nullptr) continue;
//                int i = b_timedNode->nodeID;
//                double leftT = prev_man_ptr->b_times[i];
//                for (const TimedNode *next_b_timedNode: outMapLB[b_timedNode]) {
//                    int j = next_b_timedNode->nodeID;
//                    double rightT = new_mangrove.b_times[j];
//                    TimedArc_ptr timedArc_ptr = {b_timedNode, next_b_timedNode};
//                    ttMapLB[timedArc_ptr] = G.minTT(i, j, leftT, rightT);
//                }
//            }
//        }
//    }
//    auto it_bool = mangroves[bpNode].insert(new_mangrove);
//    const Mangrove *man_ptr = &(*it_bool.first);
//    mangroveMap[bpNode][bpTime] = man_ptr;
//    return man_ptr;
//}
//
//pair<TEN::type_path, double> TEN::findLB() {
////    map<const TimedNode *, double> dp;
////    map<const TimedNode *, const TimedNode *> pred;
//
//    unordered_map<const TimedNode *, double, TimedNodeHash, TimedNodeEqual> dp;
//    unordered_map<const TimedNode *, const TimedNode *, TimedNodeHash, TimedNodeEqual> pred;
//    priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
//    pq.push({origin, 0});
//    dp[origin] = 0;
//    while (!pq.empty()) {
//        const TimedNode *ptr = pq.top().first;
//        //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpNode << ',' << ptr->bpTime << ")" << endl;
//        //cout << "DP value: " << pq.top().second << endl;
//        if (ptr == destination) break;
//        //Care double comparison
//        if (pq.top().second > dp[ptr]) {
//            pq.pop();
//            continue;
//        }
//        pq.pop();
//        for (const TimedNode *next_ptr: outMapLB[ptr]) {
//            double weight = ttMapLB[{ptr, next_ptr}];
//            auto it = dp.find(next_ptr);
//            if (it == dp.end()) dp[next_ptr] = INT_MAX;
//            if (dp[next_ptr] > dp[ptr] + weight) {
//                dp[next_ptr] = dp[ptr] + weight;
//                pq.push({next_ptr, dp[next_ptr]});
//                //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                //cout << "DP value: " << dp[next_ptr] << endl;
//                pred[next_ptr] = ptr;
//            }
//        }
//        int node = ptr->nodeID;
//        int i = ptr->bpNode;
//        if (ptr->forward) {
//            for (int j = 0; j < G.n; j++) {
//                if (bTimedNodes[node][j].empty()) continue;
//                auto next_it = bTimedNodes[node][j].upper_bound(*ptr);
//                const TimedNode *next_ptr;
//                if (j == i) {
//                    //Waiting arc to next copy
//                    if (next_it == bTimedNodes[node][j].end()) continue;
//                    else next_ptr = &(*next_it);
//                } else {
//                    //Waiting arc from i-FSPT to j-BSPT
//                    if (next_it == bTimedNodes[node][j].begin()) next_ptr = &(*next_it);
//                    else if (next_it == bTimedNodes[node][j].end()) continue;
//                    else {
//                        auto prev_it = prev(next_it, 1);
//                        if (mangroveMap[j][prev_it->bpTime]->resolved) next_ptr = &(*next_it);
//                        else next_ptr = &(*prev_it);
//                    }
//                }
//                if (next_ptr == nullptr) continue;
//                auto it = dp.find(next_ptr);
//                if (it == dp.end()) dp[next_ptr] = INT_MAX;
//                if (dp[next_ptr] > dp[ptr]) {
//                    dp[next_ptr] = dp[ptr];
//                    pq.push({next_ptr, dp[next_ptr]});
//                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                    //cout << "DP value: " << dp[next_ptr] << endl;
//                    pred[next_ptr] = ptr;
//                }
//            }
//        } else {
//            if (i == node) {
//                //Waiting arc from BSPT to FSPT through node
//                auto next_it = fTimedNodes[node][i].lower_bound(*ptr);
//                const TimedNode *next_ptr;
//                assert(next_it != fTimedNodes[node][i].end());
//                next_ptr = &(*next_it);
//                auto it = dp.find(next_ptr);
//                if (it == dp.end()) dp[next_ptr] = INT_MAX;
//                if (dp[next_ptr] > dp[ptr]) {
//                    dp[next_ptr] = dp[ptr];
//                    pq.push({next_ptr, dp[next_ptr]});
//                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                    //cout << "DP value: " << dp[next_ptr] << endl;
//                    pred[next_ptr] = ptr;
//                }
//            }
//        }
//    }
//    type_path path;
//    const TimedNode *curr = destination;
//    while (curr != origin) {
//        path.push_back(curr);
//        curr = pred[curr];
//    }
//    path.push_back(curr);
//    reverse(path.begin(), path.end());
//    return {path, dp[destination] - dp[origin]};
//}
//
//pair<TEN::type_path, double> TEN::findUB() {
//    unordered_map<const TimedNode *, double, TimedNodeHash, TimedNodeEqual> dp;
//    unordered_map<const TimedNode *, const TimedNode *, TimedNodeHash, TimedNodeEqual> pred;
//    priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
//    pq.push({origin, 0});
//    dp[origin] = 0;
//    while (!pq.empty()) {
//        const TimedNode *ptr = pq.top().first;
//        //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpNode << ',' << ptr->bpTime << ")" << endl;
//        //cout << "DP value: " << pq.top().second << endl;
//        if (ptr == destination) break;
//        //Care double comparison
//        if (pq.top().second > dp[ptr]) {
//            pq.pop();
//            continue;
//        }
//        pq.pop();
//        for (const TimedNode *next_ptr: outMapUB[ptr]) {
//            double weight = ttMapUB[{ptr, next_ptr}];
//            auto it = dp.find(next_ptr);
//            if (it == dp.end()) dp[next_ptr] = INT_MAX;
//            if (dp[next_ptr] > dp[ptr] + weight) {
//                dp[next_ptr] = dp[ptr] + weight;
//                pq.push({next_ptr, dp[next_ptr]});
//                //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                //cout << "DP value: " << dp[next_ptr] << endl;
//                pred[next_ptr] = ptr;
//            }
//        }
//        int node = ptr->nodeID;
//        int i = ptr->bpNode;
//        if (ptr->forward) {
//            for (int j = 0; j < G.n; j++) {
//                if (bTimedNodes[node][j].empty()) continue;
//                auto next_it = bTimedNodes[node][j].upper_bound(*ptr);
//                const TimedNode *next_ptr;
//                if (j == i) {
//                    //Waiting arc to next copy
//                    if (next_it == bTimedNodes[node][j].end()) continue;
//                    else next_ptr = &(*next_it);
//                } else {
//                    //Waiting arc from i-FSPT to j-BSPT
//                    if (next_it == bTimedNodes[node][j].begin()) next_ptr = &(*next_it);
//                    else if (next_it == bTimedNodes[node][j].end()) continue;
//                    next_ptr = &(*next_it);
//                }
//                if (next_ptr == nullptr) continue;
//                auto it = dp.find(next_ptr);
//                if (it == dp.end()) dp[next_ptr] = INT_MAX;
//                if (dp[next_ptr] > dp[ptr]) {
//                    dp[next_ptr] = dp[ptr];
//                    pq.push({next_ptr, dp[next_ptr]});
//                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                    //cout << "DP value: " << dp[next_ptr] << endl;
//                    pred[next_ptr] = ptr;
//                }
//            }
//        } else {
//            if (i == node) {
//                //Waiting arc from BSPT to FSPT through node
//                auto next_it = fTimedNodes[node][i].lower_bound(*ptr);
//                const TimedNode *next_ptr;
//                assert(next_it != fTimedNodes[node][i].end());
//                next_ptr = &(*next_it);
//                auto it = dp.find(next_ptr);
//                if (it == dp.end()) dp[next_ptr] = INT_MAX;
//                if (dp[next_ptr] > dp[ptr]) {
//                    dp[next_ptr] = dp[ptr];
//                    pq.push({next_ptr, dp[next_ptr]});
//                    //cout << "Pushing: (" << next_ptr->nodeID << ',' << next_ptr->nodeID << ") from: (" << next_ptr->bpNode << ',' << next_ptr->bpTime << ")" << endl;
//                    //cout << "DP value: " << dp[next_ptr] << endl;
//                    pred[next_ptr] = ptr;
//                }
//            }
//        }
//    }
//    type_path path;
//    const TimedNode *curr = destination;
//    while (curr != origin) {
//        path.push_back(curr);
//        curr = pred[curr];
//    }
//    path.push_back(curr);
//    reverse(path.begin(), path.end());
//    return {path, dp[destination] - dp[origin]};
//}
//
//vector<TypeBP> TEN::findBP(const Mangrove *curr, int addmult, int option) {
//    vector<TypeBP> mult_bps;
//    int bpNode = curr->bpNode;
//    int bpTime = -1;
//    auto next_it = mangroves[curr->bpNode].upper_bound(*curr);
//    if (next_it == mangroves[curr->bpNode].end()) return mult_bps;
//    //Check if any breakpoints in between
//    int leftT = curr->bpTime;
//    int rightT = next_it->bpTime;
//    leftT++, rightT--;
//    if (leftT <= rightT) {
//        if (addmult == 4) {
//            stack<pair<int, int>> st;
//            st.push({leftT, rightT});
//            while (!st.empty()) {
//                pair<int, int> interval = st.top();
//                st.pop();
//                leftT = interval.first;
//                rightT = interval.second;
//                if (G.outMap[bpNode].empty()) {
//                    if (option == 0) bpTime = G.min_idx_TT(G.inMap[bpNode][0], bpNode, leftT, rightT); //min
//                    else if (option == 1) bpTime = (leftT + rightT) / 2; //med
//                    else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
//                } else {
//                    if (option == 0) bpTime = G.min_idx_TT(bpNode, G.outMap[bpNode][0], leftT, rightT); //min
//                    else if (option == 1) bpTime = (leftT + rightT) / 2; //med
//                    else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
//                }
//                if (bpTime - 1 > leftT) {
//                    st.push({leftT, bpTime - 2});
//                }
//                if (bpTime + 1 < rightT) {
//                    st.push({bpTime + 2, rightT});
//                }
//                mult_bps.push_back({bpNode, bpTime});
//                if (bpTime - 1 >= leftT) mult_bps.push_back({bpNode, bpTime - 1});
//                if (bpTime + 1 <= rightT) mult_bps.push_back({bpNode, bpTime + 1});
//            }
//        } else if (addmult == 5) {
//            for (int i = leftT; i <= rightT; i++) {
//                mult_bps.push_back({bpNode, i});
//            }
//        } else {
//            if (G.outMap[bpNode].empty()) {
//                if (option == 0) bpTime = G.min_idx_TT(G.inMap[bpNode][0], bpNode, leftT, rightT); //min
//                else if (option == 1) bpTime = (leftT + rightT) / 2; //med
//                else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
//            } else {
//                if (option == 0) bpTime = G.min_idx_TT(bpNode, G.outMap[bpNode][0], leftT, rightT); //min
//                else if (option == 1) bpTime = (leftT + rightT) / 2; //med
//                else if (option == 2) bpTime = leftT + rand() % (rightT - leftT + 1); //rand
//            }
//            mult_bps.push_back({bpNode, bpTime});
//            if (addmult >= 1 && bpTime - 1 >= leftT) mult_bps.push_back({bpNode, bpTime - 1});
//            if (addmult >= 2 && bpTime + 1 <= rightT) mult_bps.push_back({bpNode, bpTime + 1});
//            if (addmult >= 3 && leftT < bpTime - 1) mult_bps.push_back({bpNode, leftT});
//        }
//    }
//    return mult_bps;
//}
//
//set<TypeBP> TEN::findBP(type_path &path) {
//    set<TypeBP> bps;
//    for (int i = 1; i < path.size(); i++) {
//        const TimedNode *timedNode = path[i];
//        if (path[i]->bpNode != path[i - 1]->bpNode) {
//            //cout << "Not Investigating: (" << timedNode->bpNode << ',' << timedNode->bpTime << ')' << endl;
//            continue;
//        }
//        //cout << "Investigating: (" << timedNode->bpNode << ',' << timedNode->bpTime << ')' << endl;
//        vector<TypeBP> mult_bp = findBP(mangroveMap[timedNode->bpNode][timedNode->bpTime]);
//        for (TypeBP bp: mult_bp) {
//            bps.insert(bp);
//        }
//    }
//    return bps;
//}
//
//bool TEN::isResolved(type_path &path) {
//    bool wait = 1;
//    for (int i = 1; i < path.size(); i++) {
//        if (path[i]->nodeID == path[i - 1]->nodeID) {
//            wait = 1;
//        } else if (wait) {
//            wait = 0;
//            //cout << path[i]->bpNode << path[i]->bpTime << endl;
//            if (!mangroveMap[path[i]->bpNode][path[i]->bpTime]->resolved) {
//                return false;
//            }
//        }
//    }
//    return true;
//}
//
//Output TEN::findMTT() {
//    int bpExplored = 0;
//    int iter = 0;
//    auto start = chrono::high_resolution_clock::now();
//    auto start2 = chrono::high_resolution_clock::now();
//    const Mangrove *lastMangrove = addMangrove(G.endN, G.endT);
//    auto stop2 = chrono::high_resolution_clock::now();
//    auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//    bpExplored++;
//    destination = lastMangrove->b_nodes[G.endN];
//    for (int i = 0; i < G.n; i++) {
//        if (i == G.endN) continue;
//        start2 = chrono::high_resolution_clock::now();
//        addMangrove(i, floor(lastMangrove->b_times[i]));
//        stop2 = chrono::high_resolution_clock::now();
//        duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//        bpExplored++;
//    }
//    start2 = chrono::high_resolution_clock::now();
//    const Mangrove *firstMangrove = addMangrove(G.startN, G.startT);
//    stop2 = chrono::high_resolution_clock::now();
//    duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//    bpExplored++;
//    origin = firstMangrove->f_nodes[G.startN];
//    for (int i = 0; i < G.n; i++) {
//        if (i == G.startN) continue;
//        start2 = chrono::high_resolution_clock::now();
//        addMangrove(i, ceil(firstMangrove->f_times[i]));
//        stop2 = chrono::high_resolution_clock::now();
//        duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//        bpExplored++;
//    }
//    //Add breakpoints associated with minimums or points near minimums
//    if (false) {
//        for (int i = 0; i < G.n; i++) {
//            int j = G.outMap[i][0];
//            int leftT = ceil(firstMangrove->f_times[i]) + 1;
//            int rightT = floor(lastMangrove->b_times[i]) - 1;
//            double reftime = G.minttMap[{i, j}][leftT][rightT];
//            queue<vector<int>> q;
//            q.push({G.min_idx_ttMap[{i, j}][leftT][rightT], leftT, rightT});
//            while (!q.empty()) {
//                vector<int> v = q.front();
//                int t = v[0];
//                leftT = v[1];
//                rightT = v[2];
//                q.pop();
//                start2 = chrono::high_resolution_clock::now();
//                addMangrove(i, t);
//                if (t > leftT) {
//                    addMangrove(i, t - 1);
//                    bpExplored++;
//                }
//                if (t - 1 > leftT) {
//                    addMangrove(i, t - 2);
//                    bpExplored++;
//                }
//                if (t < rightT) {
//                    addMangrove(i, t + 1);
//                    bpExplored++;
//                }
//                if (t + 1 < rightT) {
//                    addMangrove(i, t + 2);
//                    bpExplored++;
//                }
//                if (t + 2 < rightT) {
//                    addMangrove(i, t + 3);
//                    bpExplored++;
//                }
//                stop2 = chrono::high_resolution_clock::now();
//                duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//                bpExplored++;
//                if (leftT < t - 2 && G.minttMap[{i, j}][leftT][t - 3] < 1.03 * reftime) {
//                    q.push({G.min_idx_ttMap[{i, j}][leftT][t - 3], leftT, t - 3});
//                }
//                if (rightT > t + 3 && G.minttMap[{i, j}][t + 4][rightT] < 1.03 * reftime) {
//                    q.push({G.min_idx_ttMap[{i, j}][t + 4][rightT], t + 4, rightT});
//                }
//            }
//        }
//    }
//    //printCurrentMangroves();
//    auto start3 = chrono::high_resolution_clock::now();
//    pair<type_path, double> pathLB = findLB();
//    auto stop3 = chrono::high_resolution_clock::now();
//    auto duration3 = chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
//    iter++;
//    //pair<type_path, double> pathUB = findUB();
//    //printPathInfo(pathLB, pathUB);
//    //printCurrentMangroves();
//    while (!isResolved(pathLB.first)) {
//        set<TypeBP> next_bps = findBP(pathLB.first);
//        //printBPtoAdd(next_bps);
//        for (TypeBP next_bp: next_bps) {
//            start2 = chrono::high_resolution_clock::now();
//            addMangrove(next_bp.first, next_bp.second);
//            stop2 = chrono::high_resolution_clock::now();
//            duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
//            bpExplored++;
//        }
//        start3 = chrono::high_resolution_clock::now();
//        pathLB = findLB();
//        stop3 = chrono::high_resolution_clock::now();
//        //cout << chrono::duration_cast<chrono::milliseconds>(stop3 - start3).count() << endl;
//        duration3 += chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
//        //pair<type_path, double> pathUB = findUB();
//        //printPathInfo(pathLB, pathUB);
//        iter++;
//    }
//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
//    cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
//    cout << "Time taken by adding mangrove: " << duration2.count() << " milliseconds" << endl;
//    cout << "Time taken by calc: " << duration3.count() << " milliseconds" << endl;
//    cout << "Optimal Solution is: " << pathLB.second << endl;
//    //pair<type_path, double> pathUB = findUB();
//    //printPathInfo(pathLB, pathUB);
//    Output output;
//    output.md = 0;
//    output.ddd = 1;
//    output.bpExplored = bpExplored;
//    output.runtime = duration.count();
//    output.sp_runtime = duration3.count();
//    output.add_runtime = duration2.count();
//    output.iter = iter;
//    printLBPath(pathLB.first, output);
//    output.optVal = pathLB.second;
//    return output;
//}
//
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
