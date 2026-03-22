#include <iostream>
#include <vector>
#include <string>

using namespace std;

int n, m;
string board[55];
bool visited[55][55];
int dy[] = {0, 0, 1, -1};
int dx[] = {1, -1, 0, 0};
bool found = false;

void dfs(int r, int c, int pr, int pc, char color) {
    if (found) return;

    visited[r][c] = true;

    for (int i = 0; i < 4; i++) {
        int nr = r + dy[i];
        int nc = c + dx[i];

        if (nr < 0 || nr >= n || nc < 0 || nc >= m) continue;
        if (board[nr][nc] != color) continue;

        if (visited[nr][nc]) {
            if (nr != pr || nc != pc) {
                found = true;
                return;
            }
        } else {
            dfs(nr, nc, r, c, color);
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> n >> m;
    for (int i = 0; i < n; i++) {
        cin >> board[i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (!visited[i][j]) {
                dfs(i, j, -1, -1, board[i][j]);
                if (found) {
                    cout << "Yes" << "\n";
                    return 0;
                }
            }
        }
    }

    cout << "No" << "\n";
    return 0;
}