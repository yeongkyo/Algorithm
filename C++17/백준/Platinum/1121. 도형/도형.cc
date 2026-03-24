#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> L(n);
    for (int i = 0; i < n; i++) {
        cin >> L[i];
    }

    int k;
    cin >> k;

    sort(L.begin(), L.end());

    vector<vector<long long>> dp(k + 1, vector<long long>(50005, 0));
    dp[0][0] = 1;

    long long ans = 0;

    for (int i = 0; i < n; i++) {
        int max_len = L[i];
        for (int s = max_len + 1; s <= 50001; s++) {
            ans += dp[k - 1][s];
        }

        for (int c = k - 1; c >= 1; c--) {
            for (int s = 0; s <= 50001; s++) {
                if (dp[c - 1][s] > 0) {
                    int nxt_s = min(50001, s + L[i]);
                    dp[c][nxt_s] += dp[c - 1][s];
                }
            }
        }
    }

    cout << ans << "\n";

    return 0;
}