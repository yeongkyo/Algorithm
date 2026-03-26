#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>

using namespace std;

struct Edge {
    int u, v;
    long long s;
};

struct PQElement {
    long long target_W;
    int edge_idx;
    int side;

    bool operator>(const PQElement& other) const {
        return target_W > other.target_W;
    }
};

int parent_arr[300005];

int find_root(int x) {
    if (parent_arr[x] == x) return x;
    return parent_arr[x] = find_root(parent_arr[x]);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<long long> W(n + 1);
    for (int i = 1; i <= n; ++i) {
        cin >> W[i];
        parent_arr[i] = i;
    }

    vector<Edge> edges(m + 1);
    for (int i = 1; i <= m; ++i) {
        cin >> edges[i].u >> edges[i].v >> edges[i].s;
    }

    vector<priority_queue<PQElement, vector<PQElement>, greater<PQElement>>> pq(n + 1);
    vector<vector<long long>> target(m + 1, vector<long long>(2, -1));
    set<int> global_valid;

    auto evaluate = [&](int i) {
        int u = find_root(edges[i].u);
        int v = find_root(edges[i].v);
        if (u == v) return;

        long long sum_W = W[u] + W[v];
        if (sum_W >= edges[i].s) {
            global_valid.insert(i);
        } else {
            long long gap = edges[i].s - sum_W;
            long long add_val = (gap + 1) / 2;
            
            target[i][0] = W[u] + add_val;
            target[i][1] = W[v] + add_val;
            
            pq[u].push({target[i][0], i, 0});
            pq[v].push({target[i][1], i, 1});
        }
    };

    for (int i = 1; i <= m; ++i) {
        evaluate(i);
    }

    vector<int> ans;

    while (!global_valid.empty()) {
        int i = *global_valid.begin();
        global_valid.erase(global_valid.begin());

        int u = find_root(edges[i].u);
        int v = find_root(edges[i].v);

        if (u == v) continue;

        ans.push_back(i);

        if (pq[u].size() < pq[v].size()) swap(u, v);
        parent_arr[v] = u;
        W[u] += W[v];

        while (!pq[v].empty()) {
            pq[u].push(pq[v].top());
            pq[v].pop();
        }

        while (!pq[u].empty() && pq[u].top().target_W <= W[u]) {
            auto top = pq[u].top();
            pq[u].pop();

            if (top.target_W != target[top.edge_idx][top.side]) continue;

            evaluate(top.edge_idx);
        }
    }

    cout << ans.size() << "\n";
    for (size_t i = 0; i < ans.size(); ++i) {
        cout << ans[i] << (i + 1 == ans.size() ? "" : " ");
    }
    cout << "\n";

    return 0;
}