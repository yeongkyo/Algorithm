#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

long long dp[101];
int min_matches[] = {0, 0, 1, 7, 4, 2, 0, 8};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    for (int i = 1; i <= 100; i++) dp[i] = 1e18;
    
    dp[2] = 1; dp[3] = 7; dp[4] = 4; dp[5] = 2; dp[6] = 6; dp[7] = 8; dp[8] = 10;

    for (int i = 9; i <= 100; i++) {
        for (int j = 2; j <= 7; j++) {
            dp[i] = min(dp[i], dp[i - j] * 10 + min_matches[j]);
        }
    }

    int t;
    cin >> t;
    while (t--) {
        int n;
        cin >> n;

        cout << dp[n] << " ";

        if (n % 2 == 0) {
            for (int i = 0; i < n / 2; i++) cout << '1';
        } else {
            cout << '7';
            for (int i = 0; i < (n - 3) / 2; i++) cout << '1';
        }
        cout << "\n";
    }

    return 0;
}