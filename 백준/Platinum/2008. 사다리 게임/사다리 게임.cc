#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    int a, b, X, Y;
    cin >> a >> b >> X >> Y;

    vector<int> dp(n + 1);
    for (int j = 1; j <= n; ++j) {
        dp[j] = abs(a - j) * Y;
    }

    for (int i = 0; i < m; ++i) {
        int p;
        cin >> p;

        vector<int> next_dp = dp;
        
        next_dp[p] = min(dp[p + 1], dp[p] + X);
        next_dp[p + 1] = min(dp[p], dp[p + 1] + X);

        for (int j = 2; j <= n; ++j) {
            next_dp[j] = min(next_dp[j], next_dp[j - 1] + Y);
        }
        
        for (int j = n - 1; j >= 1; --j) {
            next_dp[j] = min(next_dp[j], next_dp[j + 1] + Y);
        }

        dp = next_dp;
    }

    cout << dp[b] << "\n";

    return 0;
}