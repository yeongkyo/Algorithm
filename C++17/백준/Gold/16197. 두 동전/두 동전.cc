#include <iostream>
#include <vector>
#include <queue>

using namespace std;

struct State {
    int y1, x1, y2, x2, cnt;
};

int n, m;
char board[20][20];
int dy[] = {-1, 1, 0, 0};
int dx[] = {0, 0, -1, 1};

bool is_out(int y, int x) {
    return (y < 0 || y >= n || x < 0 || x >= m);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> n >> m;
    int y1 = -1, x1 = -1, y2 = -1, x2 = -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> board[i][j];
            if (board[i][j] == 'o') {
                if (y1 == -1) {
                    y1 = i; x1 = j;
                } else {
                    y2 = i; x2 = j;
                }
            }
        }
    }

    queue<State> q;
    q.push({y1, x1, y2, x2, 0});

    while (!q.empty()) {
        State cur = q.front();
        q.pop();

        if (cur.cnt >= 10) continue;

        for (int i = 0; i < 4; i++) {
            int ny1 = cur.y1 + dy[i];
            int nx1 = cur.x1 + dx[i];
            int ny2 = cur.y2 + dy[i];
            int nx2 = cur.x2 + dx[i];

            bool out1 = is_out(ny1, nx1);
            bool out2 = is_out(ny2, nx2);

            if (out1 && out2) continue;
            if (out1 || out2) {
                cout << cur.cnt + 1 << "\n";
                return 0;
            }

            if (board[ny1][nx1] == '#') {
                ny1 = cur.y1;
                nx1 = cur.x1;
            }
            if (board[ny2][nx2] == '#') {
                ny2 = cur.y2;
                nx2 = cur.x2;
            }

            q.push({ny1, nx1, ny2, nx2, cur.cnt + 1});
        }
    }

    cout << -1 << "\n";
    return 0;
}