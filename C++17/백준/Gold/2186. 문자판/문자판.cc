#include <iostream>
#include <vector>
#include <string>
#include <cstring>

using namespace std;

int n, m, k;
string grid[101];
string target;
int memo[101][101][81];

int dr[] = {-1, 1, 0, 0};
int dc[] = {0, 0, -1, 1};

int solve(int r, int c, int idx) {
    if (idx == target.length() - 1) return 1;
    if (memo[r][c][idx] != -1) return memo[r][c][idx];

    int ways = 0;
    for (int i = 0; i < 4; i++) {
        for (int step = 1; step <= k; step++) {
            int nr = r + dr[i] * step;
            int nc = c + dc[i] * step;

            if (nr < 0 || nr >= n || nc < 0 || nc >= m) continue;

            if (grid[nr][nc] == target[idx + 1]) {
                ways += solve(nr, nc, idx + 1);
            }
        }
    }

    return memo[r][c][idx] = ways;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n >> m >> k)) return 0;

    for (int i = 0; i < n; i++) {
        cin >> grid[i];
    }
    cin >> target;

    memset(memo, -1, sizeof(memo));

    int total_paths = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (grid[i][j] == target[0]) {
                total_paths += solve(i, j, 0);
            }
        }
    }

    cout << total_paths << "\n";

    return 0;
}