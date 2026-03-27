#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<int> a(n), b(n), pos(n + 1);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
        pos[a[i]] = i;
    }
    for (int i = 0; i < n; i++) cin >> b[i];

    int s = pos[b[0]];
    bool ok1 = true, ok2 = true;

    for (int i = 0; i < n; i++) {
        if (a[(s + i) % n] != b[i]) ok1 = false;
        if (a[(s - i + n) % n] != b[i]) ok2 = false;
        if (!ok1 && !ok2) break;
    }

    cout << (ok1 || ok2 ? "good puzzle" : "bad puzzle") << '\n';
    return 0;
}