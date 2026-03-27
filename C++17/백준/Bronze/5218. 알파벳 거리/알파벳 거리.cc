#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;

    while (T--) {
        string a, b;
        cin >> a >> b;

        cout << "Distances:";
        for (int i = 0; i < (int)a.size(); i++) {
            int diff = b[i] - a[i];
            if (diff < 0) diff += 26;
            cout << ' ' << diff;
        }
        cout << '\n';
    }

    return 0;
}