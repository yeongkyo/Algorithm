#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int K, R;
    cin >> K >> R;

    int M;
    cin >> M;

    vector<pair<int, int>> loc(M);
    for (int i = 0; i < M; i++) cin >> loc[i].first >> loc[i].second;

    int N;
    cin >> N;

    vector<int> coverMask(N), people(N);
    long long r2 = 1LL * R * R;

    for (int i = 0; i < N; i++) {
        int x, y, s;
        cin >> x >> y >> s;
        people[i] = s;

        int mask = 0;
        for (int j = 0; j < M; j++) {
            long long dx = x - loc[j].first;
            long long dy = y - loc[j].second;
            if (dx * dx + dy * dy <= r2) mask |= (1 << j);
        }
        coverMask[i] = mask;
    }

    int ans = 0;
    int total = 1 << M;

    for (int mask = 0; mask < total; mask++) {
        if (__builtin_popcount(mask) > K) continue;

        int sum = 0;
        for (int i = 0; i < N; i++) {
            if (mask & coverMask[i]) sum += people[i];
        }
        ans = max(ans, sum);
    }

    cout << ans << '\n';
    return 0;
}