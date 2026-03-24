#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    map<string, int> name_to_id;
    int id_counter = 0;
    auto get_id = [&](const string& s) {
        if (name_to_id.find(s) == name_to_id.end()) {
            name_to_id[s] = id_counter++;
        }
        return name_to_id[s];
    };

    vector<pair<int, vector<int>>> edges_info;
    for (int i = 0; i < n; ++i) {
        string target;
        int k;
        cin >> target >> k;
        int v = get_id(target);
        vector<int> sources;
        for (int j = 0; j < k; ++j) {
            string src;
            cin >> src;
            sources.push_back(get_id(src));
        }
        edges_info.push_back({v, sources});
    }

    string query;
    cin >> query;

    int V = id_counter;
    vector<vector<int>> adj(V);
    for (const auto& info : edges_info) {
        int v = info.first;
        for (int u : info.second) {
            adj[u].push_back(v);
        }
    }

    vector<int> dfn(V, -1), low(V, -1), scc(V, -1);
    vector<int> st;
    vector<bool> in_st(V, false);
    int timer = 0, scc_cnt = 0;

    auto dfs = [&](auto& self, int u) -> void {
        dfn[u] = low[u] = timer++;
        st.push_back(u);
        in_st[u] = true;

        for (int v : adj[u]) {
            if (dfn[v] == -1) {
                self(self, v);
                low[u] = min(low[u], low[v]);
            } else if (in_st[v]) {
                low[u] = min(low[u], dfn[v]);
            }
        }

        if (low[u] == dfn[u]) {
            while (true) {
                int v = st.back();
                st.pop_back();
                in_st[v] = false;
                scc[v] = scc_cnt;
                if (u == v) break;
            }
            scc_cnt++;
        }
    };

    for (int i = 0; i < V; ++i) {
        if (dfn[i] == -1) dfs(dfs, i);
    }

    vector<int> nodes(V);
    for (int i = 0; i < V; ++i) nodes[i] = i;
    sort(nodes.begin(), nodes.end(), [&](int a, int b) {
        return scc[a] > scc[b];
    });

    vector<long long> score(V, 1);
    for (int u : nodes) {
        for (int v : adj[u]) {
            if (scc[u] != scc[v]) {
                score[v] += score[u];
            }
        }
    }

    cout << score[name_to_id[query]] << "\n";

    return 0;
}