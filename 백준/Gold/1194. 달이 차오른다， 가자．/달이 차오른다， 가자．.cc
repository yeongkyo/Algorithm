#include <bits/stdc++.h>
using namespace std;

struct State {
    int r, c, mask;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N, M;
    cin >> N >> M;
    vector<string> g(N);
    for (int i = 0; i < N; i++) cin >> g[i];

    int sr = -1, sc = -1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (g[i][j] == '0') {
                sr = i; sc = j;
                g[i][j] = '.'; // 시작칸은 빈칸 취급
            }
        }
    }

    // dist[r][c][mask] = 이동 횟수
    vector<vector<array<int, 64>>> dist(N, vector<array<int, 64>>(M));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            dist[i][j].fill(-1);

    queue<State> q;
    dist[sr][sc][0] = 0;
    q.push({sr, sc, 0});

    const int dr[4] = {-1, 1, 0, 0};
    const int dc[4] = {0, 0, -1, 1};

    while (!q.empty()) {
        auto [r, c, mask] = q.front();
        q.pop();

        if (g[r][c] == '1') {
            cout << dist[r][c][mask] << "\n";
            return 0;
        }

        for (int k = 0; k < 4; k++) {
            int nr = r + dr[k];
            int nc = c + dc[k];
            if (nr < 0 || nr >= N || nc < 0 || nc >= M) continue;

            char ch = g[nr][nc];
            if (ch == '#') continue;

            int nmask = mask;

            // 문: 열쇠 있어야 통과
            if ('A' <= ch && ch <= 'F') {
                int need = ch - 'A';
                if ((mask & (1 << need)) == 0) continue;
            }

            // 열쇠: 줍기
            if ('a' <= ch && ch <= 'f') {
                int key = ch - 'a';
                nmask |= (1 << key);
            }

            if (dist[nr][nc][nmask] != -1) continue;
            dist[nr][nc][nmask] = dist[r][c][mask] + 1;
            q.push({nr, nc, nmask});
        }
    }

    cout << -1 << "\n";
    return 0;
}
