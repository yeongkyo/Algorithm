#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string N;
    cin >> N;

    int cnt[10] = {0};
    for (char c : N) {
        cnt[c - '0']++;
    }

    int sixnine = cnt[6] + cnt[9];
    cnt[6] = (sixnine + 1) / 2;  // ceil
    cnt[9] = 0;

    int ans = 0;
    for (int i = 0; i < 10; i++) ans = max(ans, cnt[i]);

    cout << ans << "\n";
    return 0;
}
