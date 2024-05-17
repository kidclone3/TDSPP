//
// Created by delus on 3/29/24.
//

#include "Graph.h"

//Initialization
Graph::Graph(int n, int eT, int sT, int start_n, int end_n) {
    startT = sT;
    endT = eT;
    this->n = n;
    startN = start_n;
    if (end_n != -1) endN = end_n;
    else endN = n - 1;
    addNodes(n);
}

void Graph::addNode(int i) {
    Node newNode(i);
    nodes.push_back(newNode);
}

void Graph::addNodes(int n) {
    for (int i = 0; i < n; i++) {
        addNode(i);
    }
}

void Graph::addArc(int start_idx, int end_idx, vector<double> travel_times) {

    /*
    adds arc information to maps
    to get ifloor_traveltimes:
        for each integer going forwards, see where it lands, if it overshoots, assign the previous to ifloor
    to get iceil_travel_times:
        for each integer going backwards, see where it lands, if it undershoots, assign the next to iceil
    to get minttMap:


    */
    //Simple Maps
    TypeArc arc = {start_idx, end_idx};
    inMap[end_idx].push_back(start_idx);
    outMap[start_idx].push_back(end_idx);
    ttMap[arc] = travel_times;
    //ifloorMap
    int T = endT - startT + 1;
    int curr = 0;
    vector<int> ifloorTravel_times(T);
    for (int i = 0; i < T; i++) {
        while (curr <= T - 1 && curr + travel_times[curr] <= i) {
            curr++;
        }
        ifloorTravel_times[i] = curr - 1;
        //cout << "floor:";
        //cout << i << ',' << travel_times[i] << ',';
        //cout << i << ',' << ifloor_traveltimes[i] << endl;
    }
    ittfloorMap[arc] = ifloorTravel_times;
    //iceilMap
    curr = T - 1;
    vector<int> iceil_travel_times(T);
    for (int i = T - 1; i >= 0; i--) {
        while (curr >= 0 && curr + travel_times[curr] >= i) {
            curr--;
        }
        iceil_travel_times[i] = curr + 1;
        //cout << "ceil:";
        //cout << i << ',' << travel_times[i] << ',';
        //cout << i << ',' << iceil_travel_times[i] << endl;
    }
    ittceilMap[arc] = iceil_travel_times;
    //minttMap and min_idx_ttMap
    vector<vector<double>> min_tt(T, vector<double>(T, INT_MAX));
    vector<vector<int>> minIdx_tt(T, vector<int>(T, -1));
    for (int left = 0; left < T; left++) {
        double currMin = travel_times[left];
        int currIdx = left;
        for (int right = left; right < T; right++) {
            if (currMin > travel_times[right]) {
                currMin = travel_times[right];
                currIdx = right;
            }
            min_tt[left][right] = currMin;
            minIdx_tt[left][right] = currIdx;
        }
    }
    minttMap[arc] = min_tt;
    min_idx_ttMap[arc] = minIdx_tt;
}


//Functions
double Graph::TT(int start_idx, int end_idx, double start_t) { // TT = Time Travel
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

double Graph::iTT(int start_idx, int end_idx, double end_t) { // iTT = Inverse Time Travel
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

double Graph::FSP(const vector<double> &start_ts, const vector<double> &end_ts) { // Forward shortest path
    vector<double> fspt(n, endT + 1); // Forward shortest path tree
    priority_queue<type_timedNode, vector<type_timedNode>, timedNode_greater> pq;
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

vector<double> Graph::FSPT(int start_idx, double start_t, vector<TypeArc> &arcs) {
    vector<double> fspt(n, endT + 1);
    vector<int> pred(n, -1);
    priority_queue<type_timedNode, vector<type_timedNode>, timedNode_greater> pq;
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

vector<double> Graph::BSPT(int end_idx, double end_t, vector<TypeArc> &arcs) {
    vector<double> bspt(n, startT - 1); // Backward shortest path tree
    vector<int> successor(n, -1);
    priority_queue<type_timedNode, vector<type_timedNode>, timedNode_less> pq;
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

double Graph::minTT(int start_idx, int end_idx, double left_t, double right_t) {
    /* Finds departure time of minimum travel time for arc (start_idx,end_idx) in interval [left_t,right_t]*/
    left_t = min((double) endT - startT, max(0.0, left_t - startT));
    right_t = min((double) endT - startT, max(0.0, right_t - startT));
    double ans = min(TT(start_idx, end_idx, left_t), TT(start_idx, end_idx, right_t));
    return min(ans, minttMap[{start_idx, end_idx}][ceil(left_t)][floor(right_t)]);
}

int Graph::min_idx_TT(int start_idx, int end_idx, int left_t, int right_t) {
    return min_idx_ttMap[{start_idx, end_idx}][left_t][right_t];
}

bool Graph::checkFIFO() {
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