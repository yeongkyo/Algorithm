#include <bits/stdc++.h>
using namespace std;

struct Rect {
    int x1, y1, x2, y2;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    while ( (cin >> n) ) {
        if (n == 0) break;

        vector<Rect> r(n);
        for (int i = 0; i < n; i++) {
            cin >> r[i].x1 >> r[i].y1 >> r[i].x2 >> r[i].y2;
        }

        // 위상 순서가 되도록 x1 오름차순 정렬 (동률은 아무렇게나 해도 되지만 안정적으로)
        sort(r.begin(), r.end(), [](const Rect& a, const Rect& b) {
            if (a.x1 != b.x1) return a.x1 < b.x1;
            if (a.y1 != b.y1) return a.y1 < b.y1;
            if (a.x2 != b.x2) return a.x2 < b.x2;
            return a.y2 < b.y2;
        });

        vector<int> dp(n, 1);
        int ans = 1;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (r[j].x2 < r[i].x1 && r[j].y2 < r[i].y1) {
                    dp[i] = max(dp[i], dp[j] + 1);
                }
            }
            ans = max(ans, dp[i]);
        }

        cout << ans << "\n";
    }
    return 0;
}