#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const int INF = 1e9;

struct Edge {
    int to, capacity, flow, cost, rev;
};

vector<vector<Edge>> adj;

void add_edge(int from, int to, int capacity, int cost) {
    adj[from].push_back({to, capacity, 0, cost, (int)adj[to].size()});
    adj[to].push_back({from, 0, 0, -cost, (int)adj[from].size() - 1});
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    int source = 0;
    int sink = n + m + 1;
    adj.assign(sink + 1, vector<Edge>());

    vector<int> a(n + 1);
    for (int i = 1; i <= n; i++) {
        cin >> a[i];
        add_edge(m + i, sink, a[i], 0);
    }

    vector<int> b(m + 1);
    for (int i = 1; i <= m; i++) {
        cin >> b[i];
        add_edge(source, i, b[i], 0);
    }

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int cost;
            cin >> cost;
            add_edge(i, m + j, INF, cost);
        }
    }

    int total_cost = 0;

    while (true) {
        vector<int> dist(sink + 1, INF);
        vector<int> parent(sink + 1, -1);
        vector<int> edge_idx(sink + 1, -1);
        vector<bool> in_queue(sink + 1, false);
        queue<int> q;

        dist[source] = 0;
        q.push(source);
        in_queue[source] = true;

        while (!q.empty()) {
            int curr = q.front();
            q.pop();
            in_queue[curr] = false;

            for (int i = 0; i < adj[curr].size(); i++) {
                Edge& edge = adj[curr][i];
                if (edge.capacity - edge.flow > 0 && dist[edge.to] > dist[curr] + edge.cost) {
                    dist[edge.to] = dist[curr] + edge.cost;
                    parent[edge.to] = curr;
                    edge_idx[edge.to] = i;
                    if (!in_queue[edge.to]) {
                        q.push(edge.to);
                        in_queue[edge.to] = true;
                    }
                }
            }
        }

        if (parent[sink] == -1) break;

        int amount = INF;
        for (int p = sink; p != source; p = parent[p]) {
            int curr = parent[p];
            int idx = edge_idx[p];
            amount = min(amount, adj[curr][idx].capacity - adj[curr][idx].flow);
        }

        for (int p = sink; p != source; p = parent[p]) {
            int curr = parent[p];
            int idx = edge_idx[p];
            int rev_idx = adj[curr][idx].rev;
            adj[curr][idx].flow += amount;
            adj[p][rev_idx].flow -= amount;
            total_cost += amount * adj[curr][idx].cost;
        }
    }

    cout << total_cost << "\n";

    return 0;
}