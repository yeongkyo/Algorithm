#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int N, C, Q;
vector<int> adj[200005];
int in[200005], out[200005], depth[200005];
long long tree[200005];
int timer;

void dfs(int u, int p, int d) {
    in[u] = ++timer;
    depth[u] = d;
    for (int v : adj[u]) {
        if (v != p) dfs(v, u, d + 1);
    }
    out[u] = timer;
}

void update(int i, int val) {
    while (i <= N) {
        tree[i] += val;
        i += (i & -i);
    }
}

long long query(int i) {
    long long res = 0;
    while (i > 0) {
        res += tree[i];
        i -= (i & -i);
    }
    return res;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> N >> C;
    for (int i = 0; i < N - 1; i++) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    dfs(C, 0, 1);

    cin >> Q;
    while (Q--) {
        int type, a;
        cin >> type >> a;
        if (type == 1) {
            update(in[a], 1);
        } else {
            long long cnt = query(out[a]) - query(in[a] - 1);
            cout << cnt * depth[a] << "\n";
        }
    }

    return 0;
}