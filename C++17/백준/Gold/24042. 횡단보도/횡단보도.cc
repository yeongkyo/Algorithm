#include <iostream>
#include <vector>
#include <queue>

using namespace std;

const long long INF = 1e18;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    long long m;
    if (!(cin >> n >> m)) return 0;

    vector<vector<pair<int, long long>>> adj(n + 1);
    for (long long i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back({v, i});
        adj[v].push_back({u, i});
    }

    vector<long long> dist(n + 1, INF);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;

    dist[1] = 0;
    pq.push({0, 1});

    while (!pq.empty()) {
        long long curr_time = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (dist[u] < curr_time) continue;

        if (u == n) {
            cout << curr_time << "\n";
            return 0;
        }

        for (auto& edge : adj[u]) {
            int v = edge.first;
            long long cycle_idx = edge.second;
            
            long long rem = curr_time % m;
            long long wait_time = (cycle_idx - rem + m) % m;
            long long next_time = curr_time + wait_time + 1;

            if (dist[v] > next_time) {
                dist[v] = next_time;
                pq.push({next_time, v});
            }
        }
    }

    return 0;
}