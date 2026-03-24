#include <iostream>
#include <vector>

using namespace std;

vector<vector<int>> adj;
vector<int> bMatch;
vector<bool> visited;

bool dfs(int a) {
    for (int b : adj[a]) {
        if (visited[b]) continue;
        visited[b] = true;
        if (bMatch[b] == -1 || dfs(bMatch[b])) {
            bMatch[b] = a;
            return true;
        }
    }
    return false;
}

void solve() {
    int n, m;
    cin >> n >> m;
    
    adj.assign(n, vector<int>());
    bMatch.assign(n, -1);
    
    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
    }
    
    int match = 0;
    for (int i = 0; i < n; ++i) {
        visited.assign(n, false);
        if (dfs(i)) match++;
    }
    
    cout << match << "\n";
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