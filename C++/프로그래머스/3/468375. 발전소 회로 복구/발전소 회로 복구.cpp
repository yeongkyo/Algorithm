#include <vector>
#include <string>
#include <queue>
#include <cmath>
#include <algorithm>

using namespace std;

int dist2D[40][40][40][40];
int dr[4] = {-1, 1, 0, 0};
int dc[4] = {0, 0, -1, 1};

int solution(int h, vector<string> grid, vector<vector<int>> panels, vector<vector<int>> seqs) {
    int n = grid.size();
    int m = grid[0].size();
    int k = panels.size();

    for (int r1 = 0; r1 < n; ++r1) {
        
        for (int c1 = 0; c1 < m; ++c1) {
            for (int r2 = 0; r2 < n; ++r2) {
                for (int c2 = 0; c2 < m; ++c2) {
                    dist2D[r1][c1][r2][c2] = 1e9;
                }
            }
            if (grid[r1][c1] == '#') continue;

            queue<pair<int, int>> q;
            dist2D[r1][c1][r1][c1] = 0;
            q.push({r1, c1});

            while (!q.empty()) {
                auto [r, c] = q.front();
                q.pop();

                for (int d = 0; d < 4; ++d) {
                    int nr = r + dr[d];
                    int nc = c + dc[d];
                    if (nr >= 0 && nr < n && nc >= 0 && nc < m && grid[nr][nc] != '#') {
                        if (dist2D[r1][c1][nr][nc] > dist2D[r1][c1][r][c] + 1) {
                            dist2D[r1][c1][nr][nc] = dist2D[r1][c1][r][c] + 1;
                            q.push({nr, nc});
                        }
                    }
                }
            }
        }
    }

    int er = -1, ec = -1;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            if (grid[r][c] == '@') {
                er = r;
                ec = c;
            }
        }
    }

    vector<int> p_floor(k), p_row(k), p_col(k);
    for (int i = 0; i < k; ++i) {
        p_floor[i] = panels[i][0] - 1;
        p_row[i] = panels[i][1] - 1;
        p_col[i] = panels[i][2] - 1;
    }

    auto get_dist = [&](int i, int j) -> int {
        if (p_floor[i] == p_floor[j]) {
            return dist2D[p_row[i]][p_col[i]][p_row[j]][p_col[j]];
        } else {
            return dist2D[p_row[i]][p_col[i]][er][ec] 
                 + abs(p_floor[i] - p_floor[j]) 
                 + dist2D[er][ec][p_row[j]][p_col[j]];
        }
    };

    vector<int> prereq(k, 0);
    for (const auto& seq : seqs) {
        int a = seq[0] - 1;
        int b = seq[1] - 1;
        prereq[b] |= (1 << a);
    }

    vector<vector<int>> dp(1 << k, vector<int>(k, 1e9));

    if (prereq[0] == 0) {
        dp[1 << 0][0] = 0;
    } else {
        for (int v = 0; v < k; ++v) {
            if (prereq[v] == 0) {
                dp[1 << v][v] = min(dp[1 << v][v], get_dist(0, v));
            }
        }
    }

    for (int mask = 1; mask < (1 << k); ++mask) {
        for (int u = 0; u < k; ++u) {
            if (dp[mask][u] >= 1e9) continue;
            if (!(mask & (1 << u))) continue;

            for (int v = 0; v < k; ++v) {
                if (mask & (1 << v)) continue;
                if ((mask & prereq[v]) == prereq[v]) {
                    int nmask = mask | (1 << v);
                    dp[nmask][v] = min(dp[nmask][v], dp[mask][u] + get_dist(u, v));
                }
            }
        }
    }

    int ans = 1e9;
    for (int u = 0; u < k; ++u) {
        ans = min(ans, dp[(1 << k) - 1][u]);
    }

    return ans;
}