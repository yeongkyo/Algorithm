#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    int ans = 0;
    for (int i = 1; i <= N; i++) {
        int x = i;
        while (x > 0) {
            int d = x % 10;
            if (d == 3 || d == 6 || d == 9) ans++;
            x /= 10;
        }
    }

    cout << ans << '\n';
    return 0;
}