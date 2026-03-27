#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N, S;
    cin >> N >> S;

    if (N <= 2) {
        cout << 0 << '\n';
        return 0;
    }

    int M = N - 2;
    vector<bitset<1001>> dp(M + 1), ndp(M + 1);
    dp[0][0] = 1;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= M; j++) ndp[j].reset();

        for (int used = 0; used <= M; used++) {
            for (int add = 0; used + add <= M; add++) {
                int t = add * (add + 1) / 2;
                if (t <= 1000) ndp[used + add] |= (dp[used] << t);
            }
        }

        dp.swap(ndp);
    }

    cout << (dp[M][S] ? 1 : 0) << '\n';
    return 0;
}