#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

struct Edge {
    int v, id;
};

int N, M, A, B;
vector<vector<Edge>> adj;
vector<bool> is_bridge;
vector<int> dfn, low;
int timer_cnt = 0;

void dfs_bridge(int u, int p_edge) {
    dfn[u] = low[u] = ++timer_cnt;
    for (auto& edge : adj[u]) {
        int v = edge.v;
        int id = edge.id;
        if (id == p_edge) continue;

        if (dfn[v] != -1) {
            low[u] = min(low[u], dfn[v]);
        } else {
            dfs_bridge(v, id);
            low[u] = min(low[u], low[v]);
            if (low[v] > dfn[u]) {
                is_bridge[id] = true;
            }
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M >> A >> B)) return 0;

    adj.resize(N);
    is_bridge.assign(M, false);
    dfn.assign(N, -1);
    low.assign(N, -1);

    vector<pair<int, int>> edges(M);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back({v, i});
        adj[v].push_back({u, i});
        edges[i] = {u, v};
    }

    for (int i = 0; i < N; ++i) {
        if (dfn[i] == -1) {
            dfs_bridge(i, -1);
        }
    }

    vector<int> comp(N, -1);
    int num_comps = 0;
    for (int i = 0; i < N; ++i) {
        if (comp[i] == -1) {
            queue<int> q;
            q.push(i);
            comp[i] = num_comps;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (auto& edge : adj[u]) {
                    if (is_bridge[edge.id]) continue;
                    int v = edge.v;
                    if (comp[v] == -1) {
                        comp[v] = num_comps;
                        q.push(v);
                    }
                }
            }
            num_comps++;
        }
    }

    vector<vector<int>> tree_adj(num_comps);
    for (int i = 0; i < M; ++i) {
        if (is_bridge[i]) {
            int u = comp[edges[i].first];
            int v = comp[edges[i].second];
            tree_adj[u].push_back(v);
            tree_adj[v].push_back(u);
        }
    }

    int start_comp = comp[A];
    int target_comp = comp[B];

    vector<int> parent_comp(num_comps, -1);
    vector<bool> vis(num_comps, false);
    queue<int> q;

    q.push(start_comp);
    vis[start_comp] = true;

    while (!q.empty()) {
        int c = q.front();
        q.pop();
        if (c == target_comp) break;
        for (int nxt : tree_adj[c]) {
            if (!vis[nxt]) {
                vis[nxt] = true;
                parent_comp[nxt] = c;
                q.push(nxt);
            }
        }
    }

    vector<bool> on_path(num_comps, false);
    int curr = target_comp;
    while (curr != -1) {
        on_path[curr] = true;
        curr = parent_comp[curr];
    }

    vector<int> safe_nodes;
    for (int i = 0; i < N; ++i) {
        if (!on_path[comp[i]]) {
            safe_nodes.push_back(i);
        }
    }

    cout << safe_nodes.size() << "\n";
    for (int node : safe_nodes) {
        cout << node << "\n";
    }

    return 0;
}