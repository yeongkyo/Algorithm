#include <bits/stdc++.h>
using namespace std;

int round_half_up(int x, int y) {
    return (x * 2 + y) / (2 * y);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    if (n == 0) {
        cout << 0 << '\n';
        return 0;
    }

    vector<int> a(n);
    for (int i = 0; i < n; i++) cin >> a[i];

    sort(a.begin(), a.end());

    int cut = round_half_up(15 * n, 100);
    int sum = 0;
    for (int i = cut; i < n - cut; i++) sum += a[i];

    int cnt = n - 2 * cut;
    cout << round_half_up(sum, cnt) << '\n';

    return 0;
}