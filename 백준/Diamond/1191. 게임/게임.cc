#include <iostream>
#include <vector>
#include <string>
#include <queue>

using namespace std;

const int INF = 1e9;
int dr[] = {-1, 1, 0, 0};
int dc[] = {0, 0, -1, 1};

void solve() {
    int n;
    cin >> n;
    vector<string> grid(n);
    int ar = -1, ac = -1, br = -1, bc = -1;
    for (int i = 0; i < n; ++i) {
        cin >> grid[i];
        for (int j = 0; j < n; ++j) {
            if (grid[i][j] == 'A') { ar = i; ac = j; }
            if (grid[i][j] == 'B') { br = i; bc = j; }
        }
    }

    vector<vector<int>> distA(n, vector<int>(n, INF));
    vector<vector<int>> distB(n, vector<int>(n, INF));

    auto bfs = [&](int sr, int sc, vector<vector<int>>& dist) {
        queue<pair<int, int>> q;
        q.push({sr, sc});
        dist[sr][sc] = 0;
        while (!q.empty()) {
            auto [r, c] = q.front(); q.pop();
            for (int i = 0; i < 4; ++i) {
                int nr = r + dr[i];
                int nc = c + dc[i];
                if (nr >= 0 && nr < n && nc >= 0 && nc < n && grid[nr][nc] != '#') {
                    if (dist[nr][nc] == INF) {
                        dist[nr][nc] = dist[r][c] + 1;
                        q.push({nr, nc});
                    }
                }
            }
        }
    };

    bfs(ar, ac, distA);
    bfs(br, bc, distB);

    int d = distA[br][bc];
    if (d == INF || d % 2 != 0) {
        cout << "A\n";
        return;
    }

    int k = d / 2;
    vector<vector<pair<int, int>>> L_A(k + 1);
    vector<vector<pair<int, int>>> L_B(k + 1);
    vector<vector<int>> id_A(n, vector<int>(n, -1));
    vector<vector<int>> id_B(n, vector<int>(n, -1));

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (distA[r][c] + distB[r][c] == d) {
                int da = distA[r][c];
                if (da <= k) {
                    id_A[r][c] = L_A[da].size();
                    L_A[da].push_back({r, c});
                }
                int db = distB[r][c];
                if (db <= k) {
                    id_B[r][c] = L_B[db].size();
                    L_B[db].push_back({r, c});
                }
            }
        }
    }

    vector<vector<char>> next_win(L_A[k].size(), vector<char>(L_B[k].size(), 0));
    for (size_t u = 0; u < L_A[k].size(); ++u) {
        for (size_t v = 0; v < L_B[k].size(); ++v) {
            if (L_A[k][u] == L_B[k][v]) {
                next_win[u][v] = 1;
            }
        }
    }

    for (int i = k - 1; i >= 0; --i) {
        vector<vector<int>> adj_A(L_A[i].size());
        for (size_t u = 0; u < L_A[i].size(); ++u) {
            int r = L_A[i][u].first, c = L_A[i][u].second;
            for (int d_idx = 0; d_idx < 4; ++d_idx) {
                int nr = r + dr[d_idx], nc = c + dc[d_idx];
                if (nr >= 0 && nr < n && nc >= 0 && nc < n && id_A[nr][nc] != -1 && distA[nr][nc] == i + 1) {
                    adj_A[u].push_back(id_A[nr][nc]);
                }
            }
        }

        vector<vector<int>> adj_B(L_B[i].size());
        for (size_t v = 0; v < L_B[i].size(); ++v) {
            int r = L_B[i][v].first, c = L_B[i][v].second;
            for (int d_idx = 0; d_idx < 4; ++d_idx) {
                int nr = r + dr[d_idx], nc = c + dc[d_idx];
                if (nr >= 0 && nr < n && nc >= 0 && nc < n && id_B[nr][nc] != -1 && distB[nr][nc] == i + 1) {
                    adj_B[v].push_back(id_B[nr][nc]);
                }
            }
        }

        vector<vector<char>> curr_win(L_A[i].size(), vector<char>(L_B[i].size(), 0));
        for (size_t u = 0; u < L_A[i].size(); ++u) {
            for (size_t v = 0; v < L_B[i].size(); ++v) {
                bool win = true;
                for (int nu : adj_A[u]) {
                    bool found = false;
                    for (int nv : adj_B[v]) {
                        if (next_win[nu][nv]) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        win = false;
                        break;
                    }
                }
                curr_win[u][v] = win ? 1 : 0;
            }
        }
        next_win = move(curr_win);
    }

    if (next_win[0][0]) {
        cout << "B\n";
    } else {
        cout << "A\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int t;
    if (cin >> t) {
        while (t--) {
            solve();
        }
    }
    return 0;
}