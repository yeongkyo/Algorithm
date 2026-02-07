#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<long long> w(n);
    for (int i = 0; i < n; i++) cin >> w[i];

    vector<int> parent(n, -1), proto(n, -1);
    for (int i = 1; i < n; i++) cin >> parent[i] >> proto[i]; // proto: 0/1/2

    vector<long long> dp0(n, 0), dp1(n, 0);
    for (int i = 0; i < n; i++) dp1[i] = w[i];

    for (int x = n - 1; x >= 1; x--) {
        int p = parent[x];
        int t = proto[x];

        if (t == 0) { // IAmYourFriend: leaf (x adjacent only to p in current graph)
            dp1[p] += dp0[x];
            dp0[p] += max(dp0[x], dp1[x]);
        } else if (t == 2) { // WeAreYourFriends: true twin (x adjacent to p)
            long long new0 = dp0[p] + dp0[x];
            long long new1 = max(dp1[p] + dp0[x], dp0[p] + dp1[x]);
            dp0[p] = new0;
            dp1[p] = new1;
        } else { // t == 1, MyFriendsAreYourFriends: false twin (x not adjacent to p)
            long long new0 = dp0[p] + dp0[x];
            long long new1 = max({dp1[p] + dp0[x], dp0[p] + dp1[x], dp1[p] + dp1[x]});
            dp0[p] = new0;
            dp1[p] = new1;
        }
    }

    cout << max(dp0[0], dp1[0]) << "\n";
    return 0;
}
