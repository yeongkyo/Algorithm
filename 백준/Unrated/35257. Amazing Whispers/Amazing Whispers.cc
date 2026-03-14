#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M, A, B;
    if (!(cin >> N >> M >> A >> B)) return 0;

    vector<vector<int>> adj(N * M);
    for (int i = 0; i < N * (M - 1); ++i) {
        int K;
        cin >> K;
        for (int j = 0; j < K; ++j) {
            int R;
            cin >> R;
            adj[i].push_back(R);
        }
    }

    vector<vector<int>> valid_adj(N * M);

    for (int m = 0; m < M - 1; ++m) {
        vector<vector<int>> layer_adj(N);
        for (int i = 0; i < N; ++i) {
            for (int nxt : adj[m * N + i]) {
                layer_adj[i].push_back(nxt - (m + 1) * N);
            }
        }

        vector<int> match_R(N, -1);
        vector<int> match_L(N, -1);
        vector<bool> vis(N, false);

        auto dfs_match = [&](auto& self, int u) -> bool {
            for (int v : layer_adj[u]) {
                if (vis[v]) continue;
                vis[v] = true;
                if (match_R[v] == -1 || self(self, match_R[v])) {
                    match_R[v] = u;
                    match_L[u] = v;
                    return true;
                }
            }
            return false;
        };

        for (int i = 0; i < N; ++i) {
            fill(vis.begin(), vis.end(), false);
            dfs_match(dfs_match, i);
        }

        int num_nodes = 2 * N;
        vector<vector<int>> dir_adj(num_nodes);
        for (int i = 0; i < N; ++i) {
            for (int v : layer_adj[i]) {
                if (match_L[i] == v) {
                    dir_adj[N + v].push_back(i);
                } else {
                    dir_adj[i].push_back(N + v);
                }
            }
        }

        vector<int> dfn(num_nodes, -1), low(num_nodes, -1), scc(num_nodes, -1);
        int timer = 0, scc_cnt = 0;
        vector<int> st;
        vector<bool> in_st(num_nodes, false);

        auto dfs_scc = [&](auto& self, int u) -> void {
            dfn[u] = low[u] = timer++;
            st.push_back(u);
            in_st[u] = true;
            for (int v : dir_adj[u]) {
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

        for (int i = 0; i < num_nodes; ++i) {
            if (dfn[i] == -1) dfs_scc(dfs_scc, i);
        }

        for (int i = 0; i < N; ++i) {
            int u_global = m * N + i;
            for (int v : layer_adj[i]) {
                int v_global = (m + 1) * N + v;
                if (match_L[i] == v || scc[i] == scc[N + v]) {
                    valid_adj[u_global].push_back(v_global);
                }
            }
        }
    }

    vector<bool> fwd(N * M, false);
    queue<int> q;
    q.push(A);
    fwd[A] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : valid_adj[u]) {
            if (!fwd[v]) {
                fwd[v] = true;
                q.push(v);
            }
        }
    }

    vector<vector<int>> rev_valid_adj(N * M);
    for (int u = 0; u < N * M; ++u) {
        for (int v : valid_adj[u]) {
            rev_valid_adj[v].push_back(u);
        }
    }

    vector<bool> bwd(N * M, false);
    q.push(B);
    bwd[B] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : rev_valid_adj[u]) {
            if (!bwd[v]) {
                bwd[v] = true;
                q.push(v);
            }
        }
    }

    vector<int> ans;
    for (int i = 0; i < N * M; ++i) {
        if (fwd[i] && bwd[i]) {
            ans.push_back(i);
        }
    }

    cout << ans.size() << "\n";
    for (int x : ans) {
        cout << x << "\n";
    }

    return 0;
}