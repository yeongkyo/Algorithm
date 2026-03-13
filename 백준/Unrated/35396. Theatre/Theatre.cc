#include <iostream>
#include <vector>

using namespace std;

const int MOD = 1e9 + 7;

int N, M, T;
vector<pair<int, int>> edges;
long long cnt[25];

int parent_node[25], sz[25];
struct Rollback {
    int u, v, pu, pv, sz_u, sz_v;
};
vector<Rollback> history;

int find_set(int x) {
    while (x != parent_node[x]) x = parent_node[x];
    return x;
}

bool unite(int u, int v) {
    u = find_set(u);
    v = find_set(v);
    if (u == v) return false;
    if (sz[u] < sz[v]) swap(u, v);
    history.push_back({u, v, parent_node[u], parent_node[v], sz[u], sz[v]});
    parent_node[v] = u;
    sz[u] += sz[v];
    return true;
}

void rollback() {
    auto r = history.back();
    history.pop_back();
    parent_node[r.u] = r.pu;
    parent_node[r.v] = r.pv;
    sz[r.u] = r.sz_u;
    sz[r.v] = r.sz_v;
}

void dfs(int idx, int sign, int comp) {
    if (idx == M) {
        cnt[comp] = (cnt[comp] + sign + MOD) % MOD;
        return;
    }
    
    dfs(idx + 1, sign, comp);
    
    int u = edges[idx].first, v = edges[idx].second;
    bool united = unite(u, v);
    
    dfs(idx + 1, -sign, united ? comp - 1 : comp);
    
    if (united) rollback();
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M >> T)) return 0;

    vector<long long> A(T);
    for (int i = 0; i < T; ++i) {
        cin >> A[i];
    }

    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        edges.push_back({u, v});
    }

    for (int i = 0; i < N; ++i) {
        parent_node[i] = i;
        sz[i] = 1;
    }

    dfs(0, 1, N);

    for (int i = 0; i < T; ++i) {
        long long k = A[i];
        long long ans = 0;
        long long pow_k[25];
        
        pow_k[0] = 1;
        for (int c = 1; c <= N; ++c) {
            pow_k[c] = (pow_k[c - 1] * k) % MOD;
        }

        for (int c = 1; c <= N; ++c) {
            ans = (ans + cnt[c] * pow_k[c]) % MOD;
        }
        cout << ans << "\n";
    }

    return 0;
}