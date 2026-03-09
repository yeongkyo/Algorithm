#include <iostream>
#include <vector>

using namespace std;

int main() {
    int a, b, c, n;
    cin >> a >> b >> c >> n;

    vector<bool> dp(n + 1, false);
    dp[0] = true;

    int capacities[3] = {a, b, c};

    for (int i = 0; i < 3; ++i) {
        for (int j = capacities[i]; j <= n; ++j) {
            if (dp[j - capacities[i]]) {
                dp[j] = true;
            }
        }
    }

    if (dp[n]) {
        cout << 1 << endl;
    } else {
        cout << 0 << endl;
    }

    return 0;
}