#include <iostream>
#include <vector>
#include <queue>
#include <numeric>
#include <algorithm>
#include <random>

using namespace std;

int N;
vector<vector<int>> adj;

pair<int, vector<int>> stoer_wagner(vector<vector<int>> w) {
    int n = N;
    vector<int> v(n);
    iota(v.begin(), v.end(), 0);
    int best_weight = 1e9;
    vector<int> best_partition;
    vector<vector<int>> nodes(n);
    for(int i = 0; i < n; ++i) nodes[i].push_back(i);

    while (n > 1) {
        vector<int> a(n, 0), w_sum(n, 0);
        a[0] = 1;
        int prev = 0;
        for (int i = 1; i < n; ++i) w_sum[i] = w[v[0]][v[i]];

        for (int i = 1; i < n; ++i) {
            int zj = -1;
            for (int j = 1; j < n; ++j) {
                if (!a[j] && (zj == -1 || w_sum[j] > w_sum[zj])) zj = j;
            }
            a[zj] = 1;
            if (i == n - 1) {
                if (w_sum[zj] < best_weight) {
                    best_weight = w_sum[zj];
                    best_partition = nodes[v[zj]];
                }
                for (int j = 0; j < n; ++j) {
                    w[v[prev]][v[j]] += w[v[zj]][v[j]];
                    w[v[j]][v[prev]] = w[v[prev]][v[j]];
                }
                nodes[v[prev]].insert(nodes[v[prev]].end(), nodes[v[zj]].begin(), nodes[v[zj]].end());
                v.erase(v.begin() + zj);
                n--;
            } else {
                for (int j = 1; j < n; ++j) {
                    if (!a[j]) w_sum[j] += w[v[zj]][v[j]];
                }
                prev = zj;
            }
        }
    }
    return {best_weight, best_partition};
}

void solve_11() {
    auto is_conn = [&](int color) {
        vector<bool> vis(N, false);
        queue<int> q;
        q.push(0);
        vis[0] = true;
        int count = 1;
        while(!q.empty()){
            int u = q.front(); q.pop();
            for(int v = 0; v < N; ++v) {
                if (u != v && adj[u][v] == color && !vis[v]) {
                    vis[v] = true;
                    q.push(v);
                    count++;
                }
            }
        }
        return count == N;
    };

    bool c0 = is_conn(0);
    bool c1 = is_conn(1);
    if (c0 && c1) {
        cout << 0 << "\n";
        return;
    }

    int C = c0 ? 0 : 1;
    int D = 1 - C;

    vector<int> comp_id(N, -1);
    vector<vector<int>> comps;
    for(int i = 0; i < N; ++i){
        if (comp_id[i] == -1) {
            int id = comps.size();
            comps.push_back({});
            queue<int> q;
            q.push(i);
            comp_id[i] = id;
            comps.back().push_back(i);
            while(!q.empty()){
                int u = q.front(); q.pop();
                for(int v = 0; v < N; ++v){
                    if (u != v && adj[u][v] == D && comp_id[v] == -1) {
                        comp_id[v] = id;
                        comps.back().push_back(v);
                        q.push(v);
                    }
                }
            }
        }
    }

    int k = comps.size();
    if (k == 2) {
        int s0 = comps[0].size(), s1 = comps[1].size();
        if (s0 == 1 || s1 == 1) {
            int u = (s0 == 1) ? comps[0][0] : comps[1][0];
            vector<int>& big = (s0 == 1) ? comps[1] : comps[0];
            bool has_C = false;
            for(int i = 0; i < big.size(); ++i){
                for(int j = i + 1; j < big.size(); ++j){
                    if (adj[big[i]][big[j]] == C) has_C = true;
                }
            }
            if (!has_C) {
                cout << 2 << "\n";
                cout << u + 1 << " " << big[0] + 1 << "\n";
                cout << big[0] + 1 << " " << big[1] + 1 << "\n";
                return;
            }
        }
    }

    vector<pair<int, int>> valid_edges;
    for(int i = 0; i < N; ++i){
        for(int j = i + 1; j < N; ++j){
            if (adj[i][j] == C && comp_id[i] != comp_id[j]) {
                valid_edges.push_back({i, j});
            }
        }
    }

    mt19937 rng(1337);
    int attempts = 0;
    while (attempts++ < 100000) {
        shuffle(valid_edges.begin(), valid_edges.end(), rng);
        vector<int> parent(k);
        iota(parent.begin(), parent.end(), 0);
        
        auto find_set = [&](auto& self, int v) -> int {
            if (v == parent[v]) return v;
            return parent[v] = self(self, parent[v]);
        };
        
        auto unite = [&](int a, int b) {
            a = find_set(find_set, a);
            b = find_set(find_set, b);
            if (a != b) { parent[b] = a; return true; }
            return false;
        };

        vector<pair<int, int>> tree;
        for(auto& e : valid_edges) {
            if (unite(comp_id[e.first], comp_id[e.second])) {
                tree.push_back(e);
            }
        }

        if (tree.size() == k - 1) {
            for(auto& e : tree) adj[e.first][e.second] = adj[e.second][e.first] = D;
            if (is_conn(C)) {
                cout << k - 1 << "\n";
                for(auto& e : tree) {
                    cout << e.first + 1 << " " << e.second + 1 << "\n";
                }
                return;
            }
            for(auto& e : tree) adj[e.first][e.second] = adj[e.second][e.first] = C;
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N)) return 0;

    adj.assign(N, vector<int>(N));
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            cin >> adj[i][j];
        }
    }
    
    int A, B;
    cin >> A >> B;

    if (A == 0 && B == 0) {
        cout << -1 << "\n";
        return 0;
    }

    if (A == 1 && B == 1) {
        if (N == 3) {
            cout << -1 << "\n";
            return 0;
        }
        solve_11();
        return 0;
    }

    int target_color = (A == 1) ? 1 : 0;
    vector<vector<int>> w(N, vector<int>(N, 0));
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            if (i != j && adj[i][j] == target_color) {
                w[i][j] = 1;
            }
        }
    }

    pair<int, vector<int>> cut = stoer_wagner(w);
    cout << cut.first << "\n";
    
    vector<bool> in_part(N, false);
    for(int x : cut.second) in_part[x] = true;

    for(int i = 0; i < N; ++i){
        if (in_part[i]) {
            for(int j = 0; j < N; ++j){
                if (!in_part[j] && adj[i][j] == target_color) {
                    cout << i + 1 << " " << j + 1 << "\n";
                }
            }
        }
    }

    return 0;
}