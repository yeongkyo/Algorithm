#include <iostream>
#include <vector>
#include <queue>

using namespace std;

const int INF = 1e9;

struct Edge {
    int to, capacity, flow, cost, rev;
};

void add_edge(vector<vector<Edge>>& adj, int from, int to, int capacity, int cost) {
    adj[from].push_back({to, capacity, 0, cost, (int)adj[to].size()});
    adj[to].push_back({from, 0, 0, -cost, (int)adj[from].size() - 1});
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int v, e;
    while (cin >> v >> e) {
        int source = v + 1;
        int sink = v;
        int max_nodes = 2 * v + 1;
        vector<vector<Edge>> adj(max_nodes);

        for (int i = 2; i < v; i++) {
            add_edge(adj, i, i + v, 1, 0);
        }

        for (int i = 0; i < e; i++) {
            int u, w, c;
            cin >> u >> w >> c;
            int from = u + v;
            int to = w;
            add_edge(adj, from, to, 1, c);
        }

        int total_cost = 0;

        for (int step = 0; step < 2; step++) {
            vector<int> dist(max_nodes, INF);
            vector<int> parent(max_nodes, -1);
            vector<int> edge_idx(max_nodes, -1);
            vector<bool> in_queue(max_nodes, false);
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

            for (int p = sink; p != source; p = parent[p]) {
                int curr = parent[p];
                int idx = edge_idx[p];
                int rev_idx = adj[curr][idx].rev;
                adj[curr][idx].flow += 1;
                adj[p][rev_idx].flow -= 1;
                total_cost += adj[curr][idx].cost;
            }
        }

        cout << total_cost << "\n";
    }

    return 0;
}