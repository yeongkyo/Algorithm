#include <iostream>
#include <vector>

using namespace std;

int n;
int board[35][35];
long long dp[35][35][3];

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> n;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            cin >> board[i][j];
        }
    }

    dp[1][2][0] = 1;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (board[i][j] == 1) continue;

            if (j + 1 <= n && board[i][j + 1] == 0) {
                dp[i][j + 1][0] += (dp[i][j][0] + dp[i][j][2]);
            }

            if (i + 1 <= n && board[i + 1][j] == 0) {
                dp[i + 1][j][1] += (dp[i][j][1] + dp[i][j][2]);
            }

            if (i + 1 <= n && j + 1 <= n && board[i + 1][j + 1] == 0 && board[i][j + 1] == 0 && board[i + 1][j] == 0) {
                dp[i + 1][j + 1][2] += (dp[i][j][0] + dp[i][j][1] + dp[i][j][2]);
            }
        }
    }

    cout << dp[n][n][0] + dp[n][n][1] + dp[n][n][2] << "\n";

    return 0;
}