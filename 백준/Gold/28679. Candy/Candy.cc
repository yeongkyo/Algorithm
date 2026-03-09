#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, f;
    long long t;
    if (!(cin >> n >> f >> t)) return 0;

    vector<long long> a(n);
    long long total_sum = 0;
    for (int i = 0; i < n; ++i) {
        cin >> a[i];
        total_sum += a[i];
    }

    if (total_sum < t) {
        cout << "NO\n";
        return 0;
    }

    int max_k = n * (n - 1) / 2;
    vector<vector<long long>> dp(f + 1, vector<long long>(max_k + 1, -1));
    dp[0][0] = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = min(i, f - 1); j >= 0; --j) {
            for (int k = 0; k <= max_k; ++k) {
                if (dp[j][k] != -1) {
                    int next_k = k + i - j;
                    if (next_k <= max_k) {
                        dp[j + 1][next_k] = max(dp[j + 1][next_k], dp[j][k] + a[i]);
                    }
                }
            }
        }
    }

    for (int k = 0; k <= max_k; ++k) {
        if (dp[f][k] >= t) {
            cout << k << "\n";
            return 0;
        }
    }

    cout << "NO\n";
    return 0;
}