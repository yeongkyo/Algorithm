#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    long long dp[101];

    for (int i = 1; i <= 6; ++i) {
        dp[i] = i;
    }

    for (int i = 7; i <= N; ++i) {
        dp[i] = dp[i - 1] + 1;
        for (int j = 1; j <= i - 3; ++j) {
            dp[i] = max(dp[i], dp[j] * (i - j - 1));
        }
    }

    cout << dp[N] << "\n";

    return 0;
}