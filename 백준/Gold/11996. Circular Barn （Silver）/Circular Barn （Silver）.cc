#include <bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<int> c(n);
    for (int i = 0; i < n; i++) cin >> c[i];

    const ll INF = (1LL<<62);
    ll ans = INF;

    for (int s = 0; s < n; s++) {
        deque<int> q; // start indices in the traversal order [0..n-1]
        ll energy = 0;
        bool ok = true;

        for (int k = 0; k < n; k++) {
            int i = (s + k) % n;

            // cows enter at room i; in traversal index, they "start" at k
            for (int t = 0; t < c[i]; t++) q.push_back(k);

            if (q.empty()) { ok = false; break; } // can't fill this room with forward-only cows

            int start = q.front(); q.pop_front();
            ll d = (ll)k - start;
            energy += d * d;
        }

        if (ok) ans = min(ans, energy);
    }

    cout << ans << "\n";
    return 0;
}
