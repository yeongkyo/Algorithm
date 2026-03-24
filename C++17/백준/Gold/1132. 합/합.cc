#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<long long> weight(10, 0);
    vector<bool> not_zero(10, false);

    for (int i = 0; i < n; ++i) {
        string s;
        cin >> s;
        not_zero[s[0] - 'A'] = true;
        long long w = 1;
        for (int j = s.length() - 1; j >= 0; --j) {
            weight[s[j] - 'A'] += w;
            w *= 10;
        }
    }

    int zero_idx = -1;
    long long min_weight = -1;

    for (int i = 0; i < 10; ++i) {
        if (!not_zero[i]) {
            if (min_weight == -1 || weight[i] < min_weight) {
                min_weight = weight[i];
                zero_idx = i;
            }
        }
    }

    weight[zero_idx] = -1;

    sort(weight.rbegin(), weight.rend());

    long long ans = 0;
    long long mult = 9;
    for (int i = 0; i < 9; ++i) {
        if (weight[i] > 0) {
            ans += weight[i] * mult;
        }
        mult--;
    }

    cout << ans << "\n";

    return 0;
}