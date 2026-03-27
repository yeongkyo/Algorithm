#include <iostream>
#include <vector>

using namespace std;

const int MOD = 1000000007;

int parent_node[2005];

int find_root(int i) {
    if (parent_node[i] == i)
        return i;
    return parent_node[i] = find_root(parent_node[i]);
}

void unite_nodes(int u, int v) {
    int root_u = find_root(u);
    int root_v = find_root(v);
    if (root_u != root_v) {
        parent_node[root_u] = root_v;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M;
    if (!(cin >> N >> M)) return 0;

    vector<pair<int, int>> edges(M);
    for (int i = 0; i < M; ++i) {
        cin >> edges[i].first >> edges[i].second;
    }

    vector<long long> p3(M + 1);
    p3[0] = 1;
    for (int i = 1; i <= M; ++i) {
        p3[i] = (p3[i - 1] * 3) % MOD;
    }

    for (int i = 0; i < N; ++i) {
        parent_node[i] = i;
    }

    long long ans = 0;
    for (int i = M - 1; i >= 0; --i) {
        int u = edges[i].first;
        int v = edges[i].second;

        int root_u = find_root(u);
        int root_v = find_root(v);
        int root_0 = find_root(0);
        int root_N = find_root(N - 1);

        if ((root_u == root_0 && root_v == root_N) || (root_u == root_N && root_v == root_0)) {
            ans = (ans + p3[i]) % MOD;
        } else {
            unite_nodes(u, v);
        }
    }

    cout << ans << "\n";

    return 0;
}