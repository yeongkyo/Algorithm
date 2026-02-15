#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;
    while (T--) {
        int n;
        cin >> n;
        string s;
        cin >> s;

        int i = 0;
        bool ok = true;

        while (i < n) {
            if (s[i] != 'M') { ok = false; break; }
            i++;

            int cnt = 0;
            while (i + 1 < n && s[i] == 'I' && s[i+1] == 'T') {
                cnt++;
                i += 2;
            }
            if (cnt == 0) { ok = false; break; }
        }

        cout << (ok ? "YES" : "NO") << "\n";
    }
    return 0;
}
