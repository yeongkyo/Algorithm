#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<pair<long long, long long>> pts(n);
    vector<long long> x_coords(n);
    vector<long long> y_coords(n);

    for (int i = 0; i < n; ++i) {
        cin >> pts[i].first >> pts[i].second;
        x_coords[i] = pts[i].first;
        y_coords[i] = pts[i].second;
    }

    vector<long long> ans(n + 1, -1);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long cx = x_coords[i];
            long long cy = y_coords[j];

            vector<long long> dists(n);
            for (int k = 0; k < n; ++k) {
                dists[k] = abs(pts[k].first - cx) + abs(pts[k].second - cy);
            }

            sort(dists.begin(), dists.end());

            long long sum = 0;
            for (int k = 1; k <= n; ++k) {
                sum += dists[k - 1];
                if (ans[k] == -1 || sum < ans[k]) {
                    ans[k] = sum;
                }
            }
        }
    }

    for (int k = 1; k <= n; ++k) {
        cout << ans[k] << (k == n ? "" : " ");
    }
    cout << "\n";

    return 0;
}