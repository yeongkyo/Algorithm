#include <iostream>
#include <vector>
#include <string>
#include <queue>

using namespace std;

int dx[] = {-1, 1, 0, 0};
int dy[] = {0, 0, -1, 1};

void solve() {
    int N, E;
    while (cin >> N >> E && (N != 0 || E != 0)) {
        vector<string> grid(N);
        for (int i = 0; i < N; ++i) {
            cin >> grid[i];
        }

        int V = N * E;
        vector<vector<int>> adj(V);
        vector<int> left_nodes, right_nodes;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < E; ++j) {
                if (grid[i][j] == 'X') continue;
                int u = i * E + j;
                if ((i + j) % 2 == 0) {
                    left_nodes.push_back(u);
                } else {
                    right_nodes.push_back(u);
                }
                for (int d = 0; d < 4; ++d) {
                    int ni = i + dx[d];
                    int nj = j + dy[d];
                    if (ni >= 0 && ni < N && nj >= 0 && nj < E && grid[ni][nj] != 'X') {
                        int v = ni * E + nj;
                        adj[u].push_back(v);
                    }
                }
            }
        }

        vector<int> match(V, -1);
        vector<int> dist(V, -1);

        auto bfs = [&]() {
            queue<int> q;
            for (int u : left_nodes) {
                if (match[u] == -1) {
                    dist[u] = 0;
                    q.push(u);
                } else {
                    dist[u] = -1;
                }
            }
            bool has_augmenting = false;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v : adj[u]) {
                    if (match[v] == -1) {
                        has_augmenting = true;
                    } else if (dist[match[v]] == -1) {
                        dist[match[v]] = dist[u] + 1;
                        q.push(match[v]);
                    }
                }
            }
            return has_augmenting;
        };

        auto dfs = [&](auto& self, int u) -> bool {
            for (int v : adj[u]) {
                if (match[v] == -1 || (dist[match[v]] == dist[u] + 1 && self(self, match[v]))) {
                    match[v] = u;
                    match[u] = v;
                    return true;
                }
            }
            dist[u] = -1;
            return false;
        };

        while (bfs()) {
            for (int u : left_nodes) {
                if (match[u] == -1) {
                    dfs(dfs, u);
                }
            }
        }

        vector<bool> reachable(V, false);
        queue<int> qL;
        for (int u : left_nodes) {
            if (match[u] == -1) {
                reachable[u] = true;
                qL.push(u);
            }
        }
        while (!qL.empty()) {
            int u = qL.front();
            qL.pop();
            for (int v : adj[u]) {
                int u_next = match[v];
                if (u_next != -1 && !reachable[u_next]) {
                    reachable[u_next] = true;
                    qL.push(u_next);
                }
            }
        }

        queue<int> qR;
        for (int v : right_nodes) {
            if (match[v] == -1) {
                reachable[v] = true;
                qR.push(v);
            }
        }
        while (!qR.empty()) {
            int v = qR.front();
            qR.pop();
            for (int u : adj[v]) {
                int v_next = match[u];
                if (v_next != -1 && !reachable[v_next]) {
                    reachable[v_next] = true;
                    qR.push(v_next);
                }
            }
        }

        for (int i = 0; i < N; ++i) {
            string ans = "";
            for (int j = 0; j < E; ++j) {
                if (grid[i][j] == 'X') {
                    ans += 'X';
                } else {
                    int u = i * E + j;
                    if (reachable[u]) {
                        ans += 'B';
                    } else {
                        ans += 'A';
                    }
                }
            }
            cout << ans << "\n";
        }
        cout << "\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    solve();
    return 0;
}