#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const long long INF = 1e18;

int N;
vector<vector<int>> adj;
vector<long long> A;
vector<long long> in_max;
long long ans_min = INF;

void dfs1(int u, int p) {
    in_max[u] = A[u];
    for (int v : adj[u]) {
        if (v == p) continue;
        dfs1(v, u);
        in_max[u] = max(in_max[u], in_max[v] + 1);
    }
}

void dfs2(int u, int p, long long out_val) {
    ans_min = min(ans_min, max(in_max[u], out_val));
    
    long long max1 = -1, max2 = -1;
    int v1 = -1;
    for (int v : adj[u]) {
        if (v == p) continue;
        long long val = in_max[v] + 1;
        if (val > max1) {
            max2 = max1;
            max1 = val;
            v1 = v;
        } else if (val > max2) {
            max2 = val;
        }
    }
    
    for (int v : adj[u]) {
        if (v == p) continue;
        long long pass_out = max(A[u], out_val);
        if (v == v1) {
            pass_out = max(pass_out, max2);
        } else {
            pass_out = max(pass_out, max1);
        }
        dfs2(v, u, pass_out + 1);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N)) return 0;
    
    adj.resize(N + 1);
    for (int i = 0; i < N - 1; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    A.resize(N + 1);
    for (int i = 1; i <= N; ++i) {
        cin >> A[i];
    }
    
    in_max.assign(N + 1, 0);
    
    dfs1(1, 0);
    dfs2(1, 0, 0);
    
    cout << ans_min << "\n";
    
    return 0;
}