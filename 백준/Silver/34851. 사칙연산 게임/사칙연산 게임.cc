#include <iostream>
#include <vector>

using namespace std;

typedef long long ll;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    if (!(cin >> N)) return 0;

    vector<ll> a(N + 1);
    for (int i = 0; i <= N; ++i) {
        cin >> a[i];
    }

    const ll MOD = 1e9 + 7;
    ll res = 0;

    if (a[0] == 1 || a[1] == 1) {
        res = (a[0] + a[1]) % MOD;
    } else {
        res = (a[0] * a[1]) % MOD;
    }

    for (int i = 2; i <= N; ++i) {
        if (a[i] == 1) {
            res = (res + 1) % MOD;
        } else {
            res = (res * a[i]) % MOD;
        }
    }

    cout << res << endl;

    return 0;
}