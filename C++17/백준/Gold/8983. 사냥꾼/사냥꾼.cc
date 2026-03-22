#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int m, n;
    long long l;
    if (!(cin >> m >> n >> l)) return 0;

    vector<int> hunters(m);
    for (int i = 0; i < m; ++i) {
        cin >> hunters[i];
    }

    sort(hunters.begin(), hunters.end());

    int caught_count = 0;
    for (int i = 0; i < n; ++i) {
        int a, b;
        cin >> a >> b;

        if (b > l) continue;

        auto it = lower_bound(hunters.begin(), hunters.end(), a);
        bool can_catch = false;

        if (it != hunters.end()) {
            if (abs(*it - a) + b <= l) {
                can_catch = true;
            }
        }
        
        if (!can_catch && it != hunters.begin()) {
            if (abs(*(it - 1) - a) + b <= l) {
                can_catch = true;
            }
        }

        if (can_catch) {
            caught_count++;
        }
    }

    cout << caught_count << "\n";

    return 0;
}