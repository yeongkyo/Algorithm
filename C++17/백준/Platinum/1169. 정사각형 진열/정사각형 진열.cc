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

    vector<int> r(n);
    for (int i = 0; i < n; ++i) {
        cin >> r[i];
    }

    vector<int> x(n);
    x[0] = r[0];
    for (int i = 1; i < n; ++i) {
        int max_x = 0;
        for (int j = 0; j < i; ++j) {
            max_x = max(max_x, x[j] + 2 * min(r[i], r[j]));
        }
        x[i] = max_x;
    }

    vector<double> cp;
    for (int i = 0; i < n; ++i) {
        cp.push_back(x[i] - r[i]);
        cp.push_back(x[i]);
        cp.push_back(x[i] + r[i]);
        for (int j = i + 1; j < n; ++j) {
            cp.push_back((x[i] + x[j]) / 2.0 - r[i] + r[j]);
            cp.push_back((x[i] + x[j]) / 2.0 + r[i] - r[j]);
        }
    }

    sort(cp.begin(), cp.end());

    vector<bool> vis(n, false);
    for (size_t k = 0; k < cp.size() - 1; ++k) {
        if (cp[k+1] - cp[k] < 1e-7) continue;
        double mid = (cp[k] + cp[k+1]) / 2.0;
        double max_f = -1.0;
        int best_i = -1;
        
        for (int i = 0; i < n; ++i) {
            double f = 0;
            if (mid >= x[i] - r[i] && mid <= x[i] + r[i]) {
                f = 2.0 * r[i] - abs(mid - x[i]);
            }
            if (f > max_f) {
                max_f = f;
                best_i = i;
            }
        }
        
        if (best_i != -1 && max_f > 1e-7) {
            bool ok = true;
            for (int i = 0; i < n; ++i) {
                if (i == best_i) continue;
                double f = 0;
                if (mid >= x[i] - r[i] && mid <= x[i] + r[i]) {
                    f = 2.0 * r[i] - abs(mid - x[i]);
                }
                if (abs(f - max_f) < 1e-7) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                vis[best_i] = true;
            }
        }
    }

    bool first = true;
    for (int i = 0; i < n; ++i) {
        if (vis[i]) {
            if (!first) cout << " ";
            cout << (i + 1);
            first = false;
        }
    }
    cout << "\n";

    return 0;
}