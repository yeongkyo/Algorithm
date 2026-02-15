#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    long long ans = LLONG_MIN;
    long long cur = 0;

    for (int i = 0; i < n; i++) {
        long long x;
        cin >> x;
        if (i == 0) {
            cur = x;
        } else {
            cur = max(x, cur + x);
        }
        ans = max(ans, cur);
    }

    cout << ans << "\n";
    return 0;
}
