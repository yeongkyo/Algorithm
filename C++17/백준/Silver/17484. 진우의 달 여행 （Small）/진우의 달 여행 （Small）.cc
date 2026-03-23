#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int INF = 1e9;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<vector<int>> arr(n, vector<int>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> arr[i][j];
        }
    }

    vector<vector<vector<int>>> dp(n, vector<vector<int>>(m, vector<int>(3, INF)));

    for (int j = 0; j < m; j++) {
        dp[0][j][0] = arr[0][j];
        dp[0][j][1] = arr[0][j];
        dp[0][j][2] = arr[0][j];
    }

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (j > 0) {
                dp[i][j][0] = arr[i][j] + min(dp[i - 1][j - 1][1], dp[i - 1][j - 1][2]);
            }
            dp[i][j][1] = arr[i][j] + min(dp[i - 1][j][0], dp[i - 1][j][2]);
            if (j < m - 1) {
                dp[i][j][2] = arr[i][j] + min(dp[i - 1][j + 1][0], dp[i - 1][j + 1][1]);
            }
        }
    }

    int ans = INF;
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < 3; k++) {
            ans = min(ans, dp[n - 1][j][k]);
        }
    }

    cout << ans << "\n";

    return 0;
}