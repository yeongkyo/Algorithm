#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Station {
    int x, y;
};

void solve() {
    int n;
    while (cin >> n && n != 0) {
        vector<Station> stations(n);
        vector<int> ux, uy;
        for (int i = 0; i < n; ++i) {
            cin >> stations[i].x >> stations[i].y;
            ux.push_back(stations[i].x);
            uy.push_back(stations[i].y);
        }

        sort(ux.begin(), ux.end());
        ux.erase(unique(ux.begin(), ux.end()), ux.end());

        sort(uy.begin(), uy.end());
        uy.erase(unique(uy.begin(), uy.end()), uy.end());

        if (ux.size() > uy.size()) {
            for (int i = 0; i < n; ++i) {
                swap(stations[i].x, stations[i].y);
            }
            swap(ux, uy);
        }

        int k = ux.size() - 1;
        if (k < 0) k = 0;

        sort(stations.begin(), stations.end(), [](const Station& a, const Station& b) {
            if (a.y != b.y) return a.y < b.y;
            return a.x < b.x;
        });

        vector<int> left_mask(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < k; ++j) {
                if (stations[i].x > ux[j]) {
                    left_mask[i] |= (1 << j);
                }
            }
        }

        int ans = 1e9;
        vector<int> last_y(k + 1, -1);

        int max_mask = 1 << k;
        for (int mask = 0; mask < max_mask; ++mask) {
            int pop_c = __builtin_popcount(mask); 
            
            if (pop_c >= ans) continue;

            int h_walls = 0;
            int last_point = -1;
            bool valid = true;

            for (int i = 0; i <= k; ++i) last_y[i] = -1;

            for (int i = 0; i < n; ++i) {
                int col = __builtin_popcount(mask & left_mask[i]);
                int sy = stations[i].y;
                int ly = last_y[col];

                if (ly == sy) {
                    valid = false;
                    break;
                }
                
                if (ly != -1) {
                    int L = ly + 1;
                    int R = sy - 1;
                    
                    if (last_point < L) {
                        h_walls++;
                        last_point = R; 
                    }
                }
                last_y[col] = sy;
            }

            if (valid) {
                ans = min(ans, pop_c + h_walls);
            }
        }

        cout << ans << "\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    solve();
    return 0;
}