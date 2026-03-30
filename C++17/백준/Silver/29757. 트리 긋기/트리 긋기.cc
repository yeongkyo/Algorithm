#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    vector<tuple<long long, long long, int>> p;
    p.reserve(N);

    for (int i = 1; i <= N; ++i) {
        long long x, y;
        cin >> x >> y;
        p.emplace_back(x, y, i);
    }

    sort(p.begin(), p.end(), [](const auto& a, const auto& b) {
        if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
        return get<1>(a) < get<1>(b);
    });

    for (int i = 0; i + 1 < N; ++i) {
        cout << get<2>(p[i]) << ' ' << get<2>(p[i + 1]) << '\n';
    }

    return 0;
}