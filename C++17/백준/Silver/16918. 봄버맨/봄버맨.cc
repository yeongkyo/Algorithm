#include <bits/stdc++.h>
using namespace std;

vector<string> explode(const vector<string>& g) {
    int R = g.size(), C = g[0].size();
    vector<string> res(R, string(C, 'O'));
    int dx[4] = {1, -1, 0, 0};
    int dy[4] = {0, 0, 1, -1};

    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            if (g[i][j] == 'O') {
                res[i][j] = '.';
                for (int d = 0; d < 4; d++) {
                    int ni = i + dx[d], nj = j + dy[d];
                    if (0 <= ni && ni < R && 0 <= nj && nj < C) res[ni][nj] = '.';
                }
            }
        }
    }
    return res;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int R, C, N;
    cin >> R >> C >> N;

    vector<string> grid(R);
    for (int i = 0; i < R; i++) cin >> grid[i];

    if (N == 1) {
        for (auto& row : grid) cout << row << '\n';
        return 0;
    }

    if (N % 2 == 0) {
        for (int i = 0; i < R; i++) cout << string(C, 'O') << '\n';
        return 0;
    }

    vector<string> first = explode(grid);
    vector<string> second = explode(first);

    if (N % 4 == 3) {
        for (auto& row : first) cout << row << '\n';
    } else {
        for (auto& row : second) cout << row << '\n';
    }

    return 0;
}