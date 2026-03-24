#include <iostream>
#include <vector>

using namespace std;

int n, m;
vector<vector<int>> adj;
vector<int> in, out;
int timer = 0;
vector<long long> tree;

void dfs(int u) {
    in[u] = ++timer;
    for (int v : adj[u]) {
        dfs(v);
    }
    out[u] = timer;
}

void add(int i, long long w) {
    while (i <= n) {
        tree[i] += w;
        i += (i & -i);
    }
}

long long query(int i) {
    long long sum = 0;
    while (i > 0) {
        sum += tree[i];
        i -= (i & -i);
    }
    return sum;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n >> m)) return 0;

    adj.assign(n + 1, vector<int>());
    in.assign(n + 1, 0);
    out.assign(n + 1, 0);
    tree.assign(n + 2, 0);

    for (int i = 1; i <= n; i++) {
        int p;
        cin >> p;
        if (p != -1) {
            adj[p].push_back(i);
        }
    }

    dfs(1);

    for (int i = 0; i < m; i++) {
        int type;
        cin >> type;
        if (type == 1) {
            int u;
            long long w;
            cin >> u >> w;
            add(in[u], w);
            add(out[u] + 1, -w);
        } else if (type == 2) {
            int u;
            cin >> u;
            cout << query(in[u]) << "\n";
        }
    }

    return 0;
}