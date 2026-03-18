#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<int> type(n);
    vector<int> pure_color(n, -1);

    for (int i = 0; i < n; i++) {
        int cnt = 0;
        int last_color = -1;
        for (int j = 0; j < m; j++) {
            int cards;
            cin >> cards;
            if (cards > 0) {
                cnt++;
                last_color = j;
            }
        }
        if (cnt == 0) {
            type[i] = 0;
        } else if (cnt == 1) {
            type[i] = 1;
            pure_color[i] = last_color;
        } else {
            type[i] = 2;
        }
    }

    int min_moves = 1e9;

    for (int j = 0; j < n; j++) {
        int e = 0;
        int s = 0;
        vector<bool> seen(m, false);

        for (int i = 0; i < n; i++) {
            if (i == j) continue;
            
            if (type[i] == 0) {
                e++;
            } else if (type[i] == 1) {
                int c = pure_color[i];
                if (!seen[c]) {
                    seen[c] = true;
                    s++;
                }
            }
        }

        int moves = (n - 1) - e - s;
        if (moves < min_moves) {
            min_moves = moves;
        }
    }

    cout << min_moves << "\n";
    return 0;
}