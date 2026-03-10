#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

typedef long long ll;

int dp[3001][3001];
bool has_monster[3001][3001];

void solve(int case_num) {
    int R, C, K;
    cin >> R >> C >> K;

    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            has_monster[i][j] = false;
            dp[i][j] = 0;
        }
    }

    for (int i = 0; i < K; ++i) {
        int r, c;
        cin >> r >> c;
        has_monster[r][c] = true;
    }

    ll total_safe_squares = 0;

    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            if (has_monster[i][j]) {
                dp[i][j] = 0; // 몬스터가 있으면 만들 수 있는 정사각형은 0개
            } else {
                if (i == 0 || j == 0) {
                    dp[i][j] = 1;
                } else {
                    int prev_min = min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
                    dp[i][j] = prev_min + 1;
                }
                total_safe_squares += dp[i][j];
            }
        }
    }

    cout << "Case #" << case_num << ": " << total_safe_squares << "\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int T;
    cin >> T;
    for (int i = 1; i <= T; ++i) {
        solve(i);
    }

    return 0;
}