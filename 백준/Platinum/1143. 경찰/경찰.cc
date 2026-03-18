#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> cost(n);
    for (int i = 0; i < n; i++) {
        cin >> cost[i];
    }

    vector<string> adj(n);
    for (int i = 0; i < n; i++) {
        cin >> adj[i];
    }

    vector<vector<bool>> dist(n, vector<bool>(n, false));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (adj[i][j] == 'Y') dist[i][j] = true;
        }
        dist[i][i] = true;
    }

    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] && dist[k][j]) {
                    dist[i][j] = true;
                }
            }
        }
    }

    vector<int> scc(n, -1);
    int scc_cnt = 0;
    for (int i = 0; i < n; i++) {
        if (scc[i] == -1) {
            for (int j = 0; j < n; j++) {
                if (dist[i][j] && dist[j][i]) {
                    scc[j] = scc_cnt;
                }
            }
            scc_cnt++;
        }
    }

    vector<bool> in_degree_zero(scc_cnt, true);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (adj[i][j] == 'Y' && scc[i] != scc[j]) {
                in_degree_zero[scc[j]] = false;
            }
        }
    }

    vector<bool> used(n, false);
    double total_sum = 0;
    int count = 0;

    for (int c = 0; c < scc_cnt; c++) {
        if (in_degree_zero[c]) {
            int min_cost = 2000;
            int best_node = -1;
            for (int i = 0; i < n; i++) {
                if (scc[i] == c && cost[i] < min_cost) {
                    min_cost = cost[i];
                    best_node = i;
                }
            }
            total_sum += min_cost;
            count++;
            used[best_node] = true;
        }
    }

    vector<int> remaining;
    for (int i = 0; i < n; i++) {
        if (!used[i]) {
            remaining.push_back(cost[i]);
        }
    }

    sort(remaining.begin(), remaining.end());

    for (int c : remaining) {
        if (c < total_sum / count) {
            total_sum += c;
            count++;
        } else {
            break;
        }
    }

    cout << fixed << setprecision(10) << total_sum / count << "\n";

    return 0;
}