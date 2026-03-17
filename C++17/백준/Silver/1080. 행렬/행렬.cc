#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<string> a(n), b(n);
    for (int i = 0; i < n; ++i) cin >> a[i];
    for (int i = 0; i < n; ++i) cin >> b[i];

    int cnt = 0;

    if (n >= 3 && m >= 3) {
        for (int i = 0; i < n - 2; ++i) {
            for (int j = 0; j < m - 2; ++j) {
                if (a[i][j] != b[i][j]) {
                    for (int r = i; r < i + 3; ++r) {
                        for (int c = j; c < j + 3; ++c) {
                            a[r][c] = (a[r][c] == '0' ? '1' : '0');
                        }
                    }
                    cnt++;
                }
            }
        }
    }

    bool possible = true;
    for (int i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            possible = false;
            break;
        }
    }

    if (possible) cout << cnt << "\n";
    else cout << -1 << "\n";

    return 0;
}