#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <cstring>

using namespace std;

int w, h;
string board[25];
int dist[15][15];
int d_map[25][25];
int dr[] = {-1, 1, 0, 0};
int dc[] = {0, 0, -1, 1};

int bfs(pair<int, int> start, pair<int, int> end) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) d_map[i][j] = -1;
    }

    queue<pair<int, int>> q;
    q.push(start);
    d_map[start.first][start.second] = 0;

    while (!q.empty()) {
        int r = q.front().first;
        int c = q.front().second;
        q.pop();

        if (r == end.first && c == end.second) return d_map[r][c];

        for (int i = 0; i < 4; i++) {
            int nr = r + dr[i];
            int nc = c + dc[i];

            if (nr >= 0 && nr < h && nc >= 0 && nc < w && board[nr][nc] != 'x' && d_map[nr][nc] == -1) {
                d_map[nr][nc] = d_map[r][c] + 1;
                q.push({nr, nc});
            }
        }
    }
    return -1;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    while (cin >> w >> h && (w != 0 || h != 0)) {
        vector<pair<int, int>> points;
        pair<int, int> start;

        for (int i = 0; i < h; i++) {
            cin >> board[i];
            for (int j = 0; j < w; j++) {
                if (board[i][j] == 'o') start = {i, j};
                else if (board[i][j] == '*') points.push_back({i, j});
            }
        }

        points.insert(points.begin(), start);
        int p_size = points.size();
        bool possible = true;

        for (int i = 0; i < p_size; i++) {
            for (int j = i + 1; j < p_size; j++) {
                int d = bfs(points[i], points[j]);
                if (d == -1) {
                    possible = false;
                    break;
                }
                dist[i][j] = dist[j][i] = d;
            }
            if (!possible) break;
        }

        if (!possible) {
            cout << -1 << "\n";
            continue;
        }

        vector<int> p;
        for (int i = 1; i < p_size; i++) p.push_back(i);

        int min_move = 1e9;
        do {
            int current_move = 0;
            int prev = 0;
            for (int i = 0; i < p.size(); i++) {
                current_move += dist[prev][p[i]];
                prev = p[i];
            }
            min_move = min(min_move, current_move);
        } while (next_permutation(p.begin(), p.end()));

        if (min_move == 1e9) cout << 0 << "\n";
        else cout << min_move << "\n";
    }

    return 0;
}