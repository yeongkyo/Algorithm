#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

struct State {
    int r, c, dir, cost;
    bool operator>(const State& other) const {
        return cost > other.cost;
    }
};

int dr[] = {-1, 0, 1, 0};
int dc[] = {0, 1, 0, -1};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    vector<string> grid(n);
    vector<pair<int, int>> doors;
    for (int i = 0; i < n; i++) {
        cin >> grid[i];
        for (int j = 0; j < n; j++) {
            if (grid[i][j] == '#') {
                doors.push_back({i, j});
            }
        }
    }

    int startR = doors[0].first;
    int startC = doors[0].second;
    int endR = doors[1].first;
    int endC = doors[1].second;

    vector<vector<vector<int>>> dist(n, vector<vector<int>>(n, vector<int>(4, 1e9)));
    priority_queue<State, vector<State>, greater<State>> pq;

    for (int i = 0; i < 4; i++) {
        dist[startR][startC][i] = 0;
        pq.push({startR, startC, i, 0});
    }

    int ans = 1e9;

    while (!pq.empty()) {
        State curr = pq.top();
        pq.pop();

        if (curr.cost > dist[curr.r][curr.c][curr.dir]) continue;
        if (curr.r == endR && curr.c == endC) {
            ans = min(ans, curr.cost);
            continue;
        }

        int nr = curr.r + dr[curr.dir];
        int nc = curr.c + dc[curr.dir];

        if (nr < 0 || nr >= n || nc < 0 || nc >= n || grid[nr][nc] == '*') continue;

        if (dist[nr][nc][curr.dir] > curr.cost) {
            dist[nr][nc][curr.dir] = curr.cost;
            pq.push({nr, nc, curr.dir, curr.cost});
        }

        if (grid[nr][nc] == '!') {
            int nDirs[2];
            if (curr.dir % 2 == 0) {
                nDirs[0] = 1; nDirs[1] = 3;
            } else {
                nDirs[0] = 0; nDirs[1] = 2;
            }

            for (int i = 0; i < 2; i++) {
                int nd = nDirs[i];
                if (dist[nr][nc][nd] > curr.cost + 1) {
                    dist[nr][nc][nd] = curr.cost + 1;
                    pq.push({nr, nc, nd, curr.cost + 1});
                }
            }
        }
    }

    cout << ans << "\n";

    return 0;
}