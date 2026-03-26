#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n;
    int a, b, q;
    if (!(cin >> n >> a >> b >> q)) return 0;

    vector<pair<long long, char>> hints(a);
    for (int i = 0; i < a; ++i) {
        cin >> hints[i].first >> hints[i].second;
    }

    vector<long long> y(b);
    vector<long long> h(b);
    for (int i = 0; i < b; ++i) {
        cin >> y[i] >> h[i];
    }

    auto find_root = [&](long long p) {
        if (b == 0) return p;
        while (true) {
            auto it = upper_bound(y.begin(), y.end(), p);
            if (it == y.begin()) break;
            int idx = it - y.begin() - 1;
            
            if (h[idx] == 0) break;
            
            long long delta = y[idx] - h[idx];
            long long k = (p - y[idx]) / delta + 1;
            p -= k * delta;
        }
        return p;
    };

    map<long long, char> root_char;
    for (int i = 0; i < a; ++i) {
        root_char[find_root(hints[i].first)] = hints[i].second;
    }

    string result = "";
    for (int i = 0; i < q; ++i) {
        long long z;
        cin >> z;
        long long root = find_root(z);
        if (root_char.count(root)) {
            result += root_char[root];
        } else {
            result += '?';
        }
    }

    cout << result << "\n";

    return 0;
}