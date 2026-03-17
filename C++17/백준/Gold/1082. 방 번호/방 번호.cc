#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> p(n);
    int min_cost = 55;
    int min_idx = -1;
    int min_cost_non_zero = 55;
    int min_idx_non_zero = -1;

    for (int i = 0; i < n; ++i) {
        cin >> p[i];
        if (p[i] <= min_cost) {
            min_cost = p[i];
            min_idx = i;
        }
        if (i > 0 && p[i] <= min_cost_non_zero) {
            min_cost_non_zero = p[i];
            min_idx_non_zero = i;
        }
    }

    int m;
    cin >> m;

    if (n == 1 || m < min_cost_non_zero) {
        cout << "0\n";
        return 0;
    }

    string ans = "";
    ans += to_string(min_idx_non_zero);
    m -= min_cost_non_zero;

    while (m >= min_cost) {
        ans += to_string(min_idx);
        m -= min_cost;
    }

    for (int i = 0; i < ans.length(); ++i) {
        for (int j = n - 1; j > ans[i] - '0'; --j) {
            if (m + p[ans[i] - '0'] >= p[j]) {
                m += p[ans[i] - '0'] - p[j];
                ans[i] = j + '0';
                break;
            }
        }
    }

    cout << ans << "\n";
    return 0;
}