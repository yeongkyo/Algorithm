#include <iostream>
#include <vector>
#include <string>

using namespace std;

int R, C, K;
vector<string> grid;
vector<vector<bool>> visited;
int dr[] = {-1, 1, 0, 0};
int dc[] = {0, 0, -1, 1};
int ans = 0;

void dfs(int r, int c, int dist) {
    if (r == 0 && c == C - 1) {
        if (dist == K) ans++;
        return;
    }

    for (int i = 0; i < 4; i++) {
        int nr = r + dr[i];
        int nc = c + dc[i];

        if (nr >= 0 && nr < R && nc >= 0 && nc < C) {
            if (!visited[nr][nc] && grid[nr][nc] != 'T') {
                visited[nr][nc] = true;
                dfs(nr, nc, dist + 1);
                visited[nr][nc] = false;
            }
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> R >> C >> K)) return 0;

    grid.resize(R);
    for (int i = 0; i < R; i++) {
        cin >> grid[i];
    }

    visited.assign(R, vector<bool>(C, false));
    visited[R - 1][0] = true;
    dfs(R - 1, 0, 1);

    cout << ans << "\n";

    return 0;
}