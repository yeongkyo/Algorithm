#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string S, T;
    cin >> S >> T;
    int n = (int)S.size();
    int m = (int)T.size();

    // i_h(0, j) = j
    vector<int> prev_ih(m + 1), cur_ih(m + 1), cur_iv(m + 1);
    for (int j = 0; j <= m; j++) prev_ih[j] = j;

    long long ans = 0;

    for (int i = 1; i <= n; i++) {
        cur_ih[0] = 0;
        cur_iv[0] = 0; // i_v(i, 0) = 0

        for (int j = 1; j <= m; j++) {
            if (S[i - 1] == T[j - 1]) {
                // match
                cur_ih[j] = cur_iv[j - 1];
                cur_iv[j] = prev_ih[j];
            } else {
                // mismatch
                int a = prev_ih[j];
                int b = cur_iv[j - 1];
                cur_ih[j] = (a > b ? a : b);
                cur_iv[j] = (a < b ? a : b);
            }

            // contribution: (m - j + 1) * max(0, j - i_h(i, j))
            int t = j - cur_ih[j];
            if (t > 0) ans += 1LL * (m - j + 1) * t;
        }

        prev_ih.swap(cur_ih);
    }

    cout << ans << "\n";
    return 0;
}
