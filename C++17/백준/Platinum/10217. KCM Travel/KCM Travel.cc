#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Edge {
    int to, cost, time;
};

const int INF = 1e9;
int dp[105][10005];
int min_time[105];

void solve() {
    int n, m, k;
    cin >> n >> m >> k;

    vector<vector<Edge>> adj(n + 1);
    for (int i = 0; i < k; i++) {
        int u, v, c, d;
        cin >> u >> v >> c >> d;
        adj[u].push_back({v, c, d});
    }

    for (int i = 1; i <= n; i++) {
        min_time[i] = INF;
        for (int j = 0; j <= m; j++) {
            dp[i][j] = INF;
        }
    }

    dp[1][0] = 0;

    for (int c = 0; c <= m; c++) {
        for (int u = 1; u <= n; u++) {
            if (dp[u][c] == INF) continue;
            
            if (min_time[u] <= dp[u][c]) continue;
            min_time[u] = dp[u][c];

            if (u == n) continue;

            for (auto& edge : adj[u]) {
                int nc = c + edge.cost;
                if (nc > m) continue;
                
                int nt = dp[u][c] + edge.time;
                if (nt < dp[edge.to][nc]) {
                    dp[edge.to][nc] = nt;
                }
            }
        }
    }

    int ans = INF;
    for (int c = 0; c <= m; c++) {
        ans = min(ans, dp[n][c]);
    }

    if (ans == INF) cout << "Poor KCM\n";
    else cout << ans << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t;
    if (cin >> t) {
        while (t--) {
            solve();
        }
    }

    return 0;
}