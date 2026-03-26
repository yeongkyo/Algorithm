#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

const int INF = 1e9;

void solve() {
    int N, M;
    cin >> N >> M;
    
    vector<vector<int>> adj(N, vector<int>(N, 0));
    for (int i = 0; i < M; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        --u; --v;
        adj[u][v] += w;
        adj[v][u] += w;
    }

    vector<int> nodes(N);
    iota(nodes.begin(), nodes.end(), 0);
    int min_cut = INF;

    for (int i = 1; i < N; ++i) {
        vector<int> w(nodes.size(), 0);
        vector<bool> in_a(nodes.size(), false);
        int curr = -1;

        for (size_t j = 0; j < nodes.size(); ++j) {
            int next_node = -1;
            for (size_t k = 0; k < nodes.size(); ++k) {
                if (!in_a[k] && (next_node == -1 || w[k] > w[next_node])) {
                    next_node = k;
                }
            }
            
            if (j == nodes.size() - 1) {
                min_cut = min(min_cut, w[next_node]);
                for (size_t k = 0; k < nodes.size(); ++k) {
                    adj[nodes[curr]][nodes[k]] += adj[nodes[next_node]][nodes[k]];
                    adj[nodes[k]][nodes[curr]] += adj[nodes[k]][nodes[next_node]];
                }
                nodes.erase(nodes.begin() + next_node);
            } else {
                in_a[next_node] = true;
                for (size_t k = 0; k < nodes.size(); ++k) {
                    w[k] += adj[nodes[next_node]][nodes[k]];
                }
                curr = next_node;
            }
        }
    }
    cout << min_cut << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        while (T--) {
            solve();
        }
    }
    return 0;
}