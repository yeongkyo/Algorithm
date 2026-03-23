#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Coin {
    int price;
    int count;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    for (int t = 0; t < 3; t++) {
        int n;
        if (!(cin >> n)) break;

        vector<Coin> coins(n);
        int total_sum = 0;
        for (int i = 0; i < n; i++) {
            cin >> coins[i].price >> coins[i].count;
            total_sum += coins[i].price * coins[i].count;
        }

        if (total_sum % 2 != 0) {
            cout << 0 << "\n";
            continue;
        }

        int target = total_sum / 2;
        vector<bool> dp(target + 1, false);
        dp[0] = true;

        bool possible = false;
        for (int i = 0; i < n; i++) {
            int p = coins[i].price;
            int c = coins[i].count;

            for (int k = 1; c > 0; k <<= 1) {
                int num = min(k, c);
                int current_val = num * p;

                for (int j = target; j >= current_val; j--) {
                    if (dp[j - current_val]) {
                        dp[j] = true;
                    }
                }
                c -= num;
                if (dp[target]) {
                    possible = true;
                    break;
                }
            }
            if (possible) break;
        }

        if (dp[target]) cout << 1 << "\n";
        else cout << 0 << "\n";
    }

    return 0;
}