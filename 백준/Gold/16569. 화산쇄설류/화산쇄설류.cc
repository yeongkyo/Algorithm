#include <bits/stdc++.h>
using namespace std;

static const int INF = 1e9;

struct Node {
    int t, x, y;
    bool operator>(const Node& other) const {
        return t > other.t;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int M, N, V;
    cin >> M >> N >> V;

    int X, Y;
    cin >> X >> Y;
    --X; --Y;

    vector<vector<int>> h(M, vector<int>(N));
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) cin >> h[i][j];
    }

    vector<vector<int>> cover(M, vector<int>(N, INF));
    vector<vector<char>> volcano(M, vector<char>(N, 0));

    priority_queue<Node, vector<Node>, greater<Node>> pq;

    for (int i = 0; i < V; i++) {
        int x, y, t;
        cin >> x >> y >> t;
        --x; --y;
        volcano[x][y] = 1;                 // 화산 칸은 언제나 못 감

        if (t < cover[x][y]) {             // 쇄설류 도착시간(자기 칸)은 t
            cover[x][y] = t;
            pq.push({t, x, y});
        }
    }

    // 멀티소스 다익스트라: cover[x][y] = 최소 (t_i + grid 거리)
    int dx[4] = {1,-1,0,0};
    int dy[4] = {0,0,1,-1};

    while (!pq.empty()) {
        auto [t, x, y] = pq.top();
        pq.pop();
        if (t != cover[x][y]) continue;

        for (int dir = 0; dir < 4; dir++) {
            int nx = x + dx[dir], ny = y + dy[dir];
            if (nx < 0 || nx >= M || ny < 0 || ny >= N) continue;
            int nt = t + 1;
            if (nt < cover[nx][ny]) {
                cover[nx][ny] = nt;
                pq.push({nt, nx, ny});
            }
        }
    }

    // 재상이 BFS (최단시간), 단 도착 시각 < cover 여야 함. 화산칸은 금지.
    vector<vector<int>> dist(M, vector<int>(N, INF));
    queue<pair<int,int>> q;

    if (!volcano[X][Y] && 0 < cover[X][Y]) {
        dist[X][Y] = 0;
        q.push({X, Y});
    }

    while (!q.empty()) {
        auto [x, y] = q.front();
        q.pop();

        int curT = dist[x][y];

        for (int dir = 0; dir < 4; dir++) {
            int nx = x + dx[dir], ny = y + dy[dir];
            if (nx < 0 || nx >= M || ny < 0 || ny >= N) continue;
            if (volcano[nx][ny]) continue;

            int nt = curT + 1;
            if (nt >= cover[nx][ny]) continue; // 도착할 때 이미 덮였거나 같은 시각이면 불가
            if (nt < dist[nx][ny]) {
                dist[nx][ny] = nt;
                q.push({nx, ny});
            }
        }
    }

    int bestH = -1;
    int bestT = INF;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (dist[i][j] == INF) continue;
            int height = h[i][j];
            int time = dist[i][j];

            if (height > bestH) {
                bestH = height;
                bestT = time;
            } else if (height == bestH && time < bestT) {
                bestT = time;
            }
        }
    }

    cout << bestH << " " << bestT << "\n";
    return 0;
}