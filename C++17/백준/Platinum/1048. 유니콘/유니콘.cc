#include <iostream>
#include <vector>
#include <string>

using namespace std;

const int MOD = 1000000007;
int dp[55][305][305];
int psum[305][305];

int get_sum(int r1, int c1, int r2, int c2) {
    if (r1 > r2 || c1 > c2) return 0;
    long long res = psum[r2][c2];
    res = (res - psum[r1 - 1][c2] + MOD) % MOD;
    res = (res - psum[r2][c1 - 1] + MOD) % MOD;
    res = (res + psum[r1 - 1][c1 - 1]) % MOD;
    return res;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M, L;
    if (!(cin >> N >> M >> L)) return 0;

    string word;
    cin >> word;

    vector<string> board(N);
    for (int i = 0; i < N; i++) {
        cin >> board[i];
    }

    int K = word.length();

    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= M; j++) {
            if (board[i - 1][j - 1] == word[0]) {
                dp[0][i][j] = 1;
            }
        }
    }

    for (int k = 0; k < K - 1; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= M; j++) {
                long long val = psum[i - 1][j] + psum[i][j - 1];
                val = (val - psum[i - 1][j - 1] + MOD) % MOD;
                val = (val + dp[k][i][j]) % MOD;
                psum[i][j] = val;
            }
        }

        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= M; j++) {
                if (board[i - 1][j - 1] != word[k + 1]) continue;

                long long total = 0;

                total = (total + get_sum(1, 1, i - 2, j - 2)) % MOD;
                total = (total + get_sum(1, j + 2, i - 2, M)) % MOD;
                total = (total + get_sum(i + 2, 1, N, j - 2)) % MOD;
                total = (total + get_sum(i + 2, j + 2, N, M)) % MOD;

                int dr[] = {-2, -2, 2, 2};
                int dc[] = {-2, 2, -2, 2};

                for (int d = 0; d < 4; d++) {
                    int r = i + dr[d];
                    int c = j + dc[d];
                    if (r >= 1 && r <= N && c >= 1 && c <= M) {
                        total = (total - dp[k][r][c] + MOD) % MOD;
                    }
                }

                dp[k + 1][i][j] = total;
            }
        }
    }

    long long ans = 0;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= M; j++) {
            ans = (ans + dp[K - 1][i][j]) % MOD;
        }
    }

    cout << ans << "\n";
    return 0;
}