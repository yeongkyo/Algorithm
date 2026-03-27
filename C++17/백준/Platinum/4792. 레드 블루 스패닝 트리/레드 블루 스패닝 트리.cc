#include <bits/stdc++.h>
using namespace std;

struct DSU {
    vector<int> p, r;
    DSU(int n) : p(n + 1), r(n + 1, 0) {
        for (int i = 1; i <= n; i++) p[i] = i;
    }
    int find(int x) {
        if (p[x] == x) return x;
        return p[x] = find(p[x]);
    }
    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (r[a] < r[b]) swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
        return true;
    }
};

struct Edge {
    char c;
    int u, v;
};

int countBlue(int n, const vector<Edge>& edges, bool maximizeBlue) {
    DSU dsu(n);
    int blue = 0, used = 0;

    if (maximizeBlue) {
        for (const auto& e : edges) {
            if (e.c == 'B' && dsu.unite(e.u, e.v)) {
                blue++;
                used++;
            }
        }
        for (const auto& e : edges) {
            if (e.c == 'R' && dsu.unite(e.u, e.v)) {
                used++;
            }
        }
    } else {
        for (const auto& e : edges) {
            if (e.c == 'R' && dsu.unite(e.u, e.v)) {
                used++;
            }
        }
        for (const auto& e : edges) {
            if (e.c == 'B' && dsu.unite(e.u, e.v)) {
                blue++;
                used++;
            }
        }
    }

    if (used != n - 1) return -1;
    return blue;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    while (true) {
        int n, m, k;
        cin >> n >> m >> k;
        if (n == 0 && m == 0 && k == 0) break;

        vector<Edge> edges(m);
        for (int i = 0; i < m; i++) {
            cin >> edges[i].c >> edges[i].u >> edges[i].v;
        }

        int minBlue = countBlue(n, edges, false);
        int maxBlue = countBlue(n, edges, true);

        if (minBlue != -1 && maxBlue != -1 && minBlue <= k && k <= maxBlue) {
            cout << 1 << '\n';
        } else {
            cout << 0 << '\n';
        }
    }

    return 0;
}