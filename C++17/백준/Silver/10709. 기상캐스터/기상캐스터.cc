#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int H, W;
    cin >> H >> W;

    for (int i = 0; i < H; i++) {
        string s;
        cin >> s;

        int last = -1;
        for (int j = 0; j < W; j++) {
            if (s[j] == 'c') {
                last = j;
                cout << 0;
            } else {
                if (last == -1) cout << -1;
                else cout << j - last;
            }
            if (j + 1 < W) cout << ' ';
        }
        cout << '\n';
    }

    return 0;
}