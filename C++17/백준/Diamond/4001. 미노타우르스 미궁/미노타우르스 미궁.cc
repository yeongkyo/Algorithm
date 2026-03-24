#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

struct Point {
    int x, y;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int W, H;
    if (!(cin >> W >> H)) return 0;

    vector<string> grid(H);
    for (int i = 0; i < H; i++) cin >> grid[i];

    vector<vector<int>> dp(W + 5, vector<int>(H + 5, 0));
    for (int x = W; x >= 1; x--) {
        for (int y = H; y >= 1; y--) {
            if (grid[y - 1][x - 1] == '.') {
                dp[x][y] = min({dp[x + 1][y], dp[x][y + 1], dp[x + 1][y + 1]}) + 1;
            } else {
                dp[x][y] = 0;
            }
        }
    }

    vector<vector<bool>> is_wall(W + 5, vector<bool>(H + 5, false));
    vector<vector<bool>> in_G1(W + 5, vector<bool>(H + 5, false));
    vector<vector<bool>> in_G2(W + 5, vector<bool>(H + 5, false));

    for (int x = 1; x <= W; x++) {
        for (int y = 1; y <= H; y++) {
            if (grid[y - 1][x - 1] == '#') {
                is_wall[x][y] = true;
            }
        }
    }

    queue<Point> q1, q2;

    for (int x = 2; x <= W; x++) {
        in_G1[x][0] = true;
        is_wall[x][0] = true;
        q1.push({x, 0});
    }
    for (int y = 1; y <= H - 1; y++) {
        in_G1[W + 1][y] = true;
        is_wall[W + 1][y] = true;
        q1.push({W + 1, y});
    }

    for (int x = 1; x <= W - 1; x++) {
        in_G2[x][H + 1] = true;
        is_wall[x][H + 1] = true;
        q2.push({x, H + 1});
    }
    for (int y = 2; y <= H; y++) {
        in_G2[0][y] = true;
        is_wall[0][y] = true;
        q2.push({0, y});
    }

    int dx[] = {-1, -1, -1, 0, 1, 1, 1, 0};
    int dy[] = {-1, 0, 1, 1, 1, 0, -1, -1};

    while (!q1.empty()) {
        Point u = q1.front();
        q1.pop();
        for (int i = 0; i < 8; i++) {
            int nx = u.x + dx[i];
            int ny = u.y + dy[i];
            if (nx >= 0 && nx <= W + 1 && ny >= 0 && ny <= H + 1) {
                if (is_wall[nx][ny] && !in_G1[nx][ny]) {
                    in_G1[nx][ny] = true;
                    q1.push({nx, ny});
                }
            }
        }
    }

    while (!q2.empty()) {
        Point u = q2.front();
        q2.pop();
        for (int i = 0; i < 8; i++) {
            int nx = u.x + dx[i];
            int ny = u.y + dy[i];
            if (nx >= 0 && nx <= W + 1 && ny >= 0 && ny <= H + 1) {
                if (is_wall[nx][ny] && !in_G2[nx][ny]) {
                    in_G2[nx][ny] = true;
                    q2.push({nx, ny});
                }
            }
        }
    }

    vector<vector<bool>> touch_G1(W + 5, vector<bool>(H + 5, false));
    vector<vector<bool>> touch_G2(W + 5, vector<bool>(H + 5, false));

    for (int x = 1; x <= W; x++) {
        for (int y = 1; y <= H; y++) {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    int nx = x + i;
                    int ny = y + j;
                    if (nx >= 0 && nx <= W + 1 && ny >= 0 && ny <= H + 1) {
                        if (in_G1[nx][ny]) touch_G1[x][y] = true;
                        if (in_G2[nx][ny]) touch_G2[x][y] = true;
                    }
                }
            }
        }
    }

    vector<vector<int>> pref1(W + 5, vector<int>(H + 5, 0));
    vector<vector<int>> pref2(W + 5, vector<int>(H + 5, 0));

    for (int x = 1; x <= W; x++) {
        for (int y = 1; y <= H; y++) {
            pref1[x][y] = pref1[x - 1][y] + pref1[x][y - 1] - pref1[x - 1][y - 1] + (touch_G1[x][y] ? 1 : 0);
            pref2[x][y] = pref2[x - 1][y] + pref2[x][y - 1] - pref2[x - 1][y - 1] + (touch_G2[x][y] ? 1 : 0);
        }
    }

    auto get_sum1 = [&](int x1, int y1, int x2, int y2) {
        return pref1[x2][y2] - pref1[x1 - 1][y2] - pref1[x2][y1 - 1] + pref1[x1 - 1][y1 - 1];
    };
    auto get_sum2 = [&](int x1, int y1, int x2, int y2) {
        return pref2[x2][y2] - pref2[x1 - 1][y2] - pref2[x2][y1 - 1] + pref2[x1 - 1][y1 - 1];
    };

    int min_L = 1e9;
    int best_x = -1;
    int best_y = -1;

    for (int x = 1; x <= W; x++) {
        for (int y = 1; y <= H; y++) {
            if (x == 1 && y == 1) continue;
            
            int M = dp[x][y];
            if (M == 0) continue;
            
            int L_bad = max(W - x + 1, H - y + 1);
            int MaxL = min(M, L_bad - 1);
            if (MaxL < 1) continue;

            int low = 1, high = MaxL, ans_L = -1;
            while (low <= high) {
                int mid = low + (high - low) / 2;
                int s1 = get_sum1(x, y, x + mid - 1, y + mid - 1);
                int s2 = get_sum2(x, y, x + mid - 1, y + mid - 1);
                
                if (s1 > 0 && s2 > 0) {
                    ans_L = mid;
                    high = mid - 1;
                } else {
                    low = mid + 1;
                }
            }

            if (ans_L != -1 && ans_L < min_L) {
                min_L = ans_L;
                best_x = x;
                best_y = y;
            }
        }
    }

    if (min_L > 1e8) {
        cout << "Impossible\n";
    } else {
        cout << min_L << " " << best_x << " " << best_y << "\n";
    }

    return 0;
}