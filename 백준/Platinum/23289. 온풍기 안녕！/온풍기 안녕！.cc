#include <bits/stdc++.h>
using namespace std;

struct Heater {
    int r, c, d;
    vector<tuple<int,int,int>> eff;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int R, C, K;
    cin >> R >> C >> K;

    vector<vector<int>> a(R+1, vector<int>(C+1, 0));
    vector<Heater> heaters;
    vector<pair<int,int>> checks;

    for (int i = 1; i <= R; i++) {
        for (int j = 1; j <= C; j++) {
            cin >> a[i][j];
            if (1 <= a[i][j] && a[i][j] <= 4) heaters.push_back({i, j, a[i][j], {}});
            else if (a[i][j] == 5) checks.push_back({i, j});
        }
    }

    int W;
    cin >> W;
    vector<vector<bool>> wallUp(R+1, vector<bool>(C+1, false));
    vector<vector<bool>> wallRight(R+1, vector<bool>(C+1, false));

    for (int i = 0; i < W; i++) {
        int x, y, t;
        cin >> x >> y >> t;
        if (t == 0) wallUp[x][y] = true;
        else wallRight[x][y] = true;
    }

    auto inb = [&](int r, int c) {
        return 1 <= r && r <= R && 1 <= c && c <= C;
    };

    auto noWall = [&](int r1, int c1, int r2, int c2) {
        if (r1 == r2 && c2 == c1 + 1) return !wallRight[r1][c1];
        if (r1 == r2 && c2 == c1 - 1) return !wallRight[r1][c2];
        if (c1 == c2 && r2 == r1 - 1) return !wallUp[r1][c1];
        if (c1 == c2 && r2 == r1 + 1) return !wallUp[r2][c2];
        return false;
    };

    auto buildEff = [&](Heater &h) {
        int dr = 0, dc = 0;
        if (h.d == 1) dc = 1;
        else if (h.d == 2) dc = -1;
        else if (h.d == 3) dr = -1;
        else dr = 1;

        int sr = h.r + dr, sc = h.c + dc;
        if (!inb(sr, sc)) return;

        vector<vector<bool>> vis(R+1, vector<bool>(C+1, false));
        queue<tuple<int,int,int>> q;
        vis[sr][sc] = true;
        q.push({sr, sc, 5});

        while (!q.empty()) {
            auto [x, y, v] = q.front();
            q.pop();
            h.eff.push_back({x, y, v});
            if (v == 1) continue;

            auto pushIf = [&](int nx, int ny, bool ok) {
                if (!ok || !inb(nx, ny) || vis[nx][ny]) return;
                vis[nx][ny] = true;
                q.push({nx, ny, v-1});
            };

            if (h.d == 1) {
                pushIf(x-1, y+1, inb(x-1,y) && noWall(x,y,x-1,y) && inb(x-1,y+1) && noWall(x-1,y,x-1,y+1));
                pushIf(x,   y+1, inb(x,y+1) && noWall(x,y,x,y+1));
                pushIf(x+1, y+1, inb(x+1,y) && noWall(x,y,x+1,y) && inb(x+1,y+1) && noWall(x+1,y,x+1,y+1));
            } else if (h.d == 2) {
                pushIf(x-1, y-1, inb(x-1,y) && noWall(x,y,x-1,y) && inb(x-1,y-1) && noWall(x-1,y,x-1,y-1));
                pushIf(x,   y-1, inb(x,y-1) && noWall(x,y,x,y-1));
                pushIf(x+1, y-1, inb(x+1,y) && noWall(x,y,x+1,y) && inb(x+1,y-1) && noWall(x+1,y,x+1,y-1));
            } else if (h.d == 3) {
                pushIf(x-1, y-1, inb(x,y-1) && noWall(x,y,x,y-1) && inb(x-1,y-1) && noWall(x,y-1,x-1,y-1));
                pushIf(x-1, y,   inb(x-1,y) && noWall(x,y,x-1,y));
                pushIf(x-1, y+1, inb(x,y+1) && noWall(x,y,x,y+1) && inb(x-1,y+1) && noWall(x,y+1,x-1,y+1));
            } else {
                pushIf(x+1, y-1, inb(x,y-1) && noWall(x,y,x,y-1) && inb(x+1,y-1) && noWall(x,y-1,x+1,y-1));
                pushIf(x+1, y,   inb(x+1,y) && noWall(x,y,x+1,y));
                pushIf(x+1, y+1, inb(x,y+1) && noWall(x,y,x,y+1) && inb(x+1,y+1) && noWall(x,y+1,x+1,y+1));
            }
        }
    };

    for (auto &h : heaters) buildEff(h);

    vector<vector<int>> temp(R+1, vector<int>(C+1, 0));

    auto allOk = [&]() {
        for (auto &p : checks) {
            if (temp[p.first][p.second] < K) return false;
        }
        return true;
    };

    int choco = 0;
    while (true) {
        if (choco > 100) {
            cout << 101 << "\n";
            return 0;
        }
        if (allOk()) {
            cout << choco << "\n";
            return 0;
        }

        for (auto &h : heaters) {
            for (auto &[x,y,v] : h.eff) temp[x][y] += v;
        }

        vector<vector<int>> delta(R+1, vector<int>(C+1, 0));
        for (int i = 1; i <= R; i++) {
            for (int j = 1; j <= C; j++) {
                if (j < C && !wallRight[i][j]) {
                    int a = temp[i][j], b = temp[i][j+1];
                    int d = abs(a-b)/4;
                    if (d) {
                        if (a > b) { delta[i][j] -= d; delta[i][j+1] += d; }
                        else { delta[i][j] += d; delta[i][j+1] -= d; }
                    }
                }
                if (i < R && !wallUp[i+1][j]) {
                    int a = temp[i][j], b = temp[i+1][j];
                    int d = abs(a-b)/4;
                    if (d) {
                        if (a > b) { delta[i][j] -= d; delta[i+1][j] += d; }
                        else { delta[i][j] += d; delta[i+1][j] -= d; }
                    }
                }
            }
        }
        for (int i = 1; i <= R; i++)
            for (int j = 1; j <= C; j++)
                temp[i][j] += delta[i][j];

        for (int i = 1; i <= R; i++) {
            for (int j = 1; j <= C; j++) {
                if (i == 1 || i == R || j == 1 || j == C) {
                    if (temp[i][j] > 0) temp[i][j]--;
                }
            }
        }

        choco++;
    }
}