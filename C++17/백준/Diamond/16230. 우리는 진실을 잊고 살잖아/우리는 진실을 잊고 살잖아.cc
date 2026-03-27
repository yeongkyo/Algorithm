#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

int adj[505][505];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n, m;
    if (!(cin >> n >> m)) return 0;
    
    vector<vector<int>> graph(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        u--; v--;
        graph[u].push_back(v);
        graph[v].push_back(u);
        adj[u][v] = 1;
        adj[v][u] = 1;
    }
    
    vector<int> comp_size;
    vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            int sz = 0;
            queue<int> q;
            q.push(i);
            visited[i] = true;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                sz++;
                for (int v : graph[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        q.push(v);
                    }
                }
            }
            comp_size.push_back(sz);
        }
    }
    
    long long total_pairs = (long long)n * (n - 1) / 2;
    
    if (comp_size.size() > 1) {
        int min_comp = n + 1;
        for (int sz : comp_size) {
            min_comp = min(min_comp, sz);
        }
        long long min_reveals = (long long)min_comp * (n - min_comp);
        long long max_reveals = total_pairs - comp_size.size() + 2;
        cout << min_reveals << "\n" << max_reveals << "\n";
    } else {
        int min_c = 1e9;
        vector<int> v(n);
        for (int i = 0; i < n; i++) v[i] = i;
        
        int sz = n;
        while (sz > 1) {
            vector<int> w(n, 0);
            vector<bool> vis(n, false);
            int prev = 0, last = 0;
            
            for (int i = 0; i < sz; i++) {
                int nxt = -1, max_w = -1;
                for (int j = 0; j < sz; j++) {
                    if (!vis[j] && w[j] > max_w) {
                        max_w = w[j];
                        nxt = j;
                    }
                }
                if (nxt == -1) break;
                
                last = prev;
                prev = nxt;
                vis[nxt] = true;
                
                int v_prev = v[nxt];
                for (int j = 0; j < sz; j++) {
                    if (!vis[j]) w[j] += adj[v_prev][v[j]];
                }
            }
            
            min_c = min(min_c, w[prev]);
            
            for (int j = 0; j < sz; j++) {
                if (j == last || j == prev) continue;
                adj[v[last]][v[j]] += adj[v[prev]][v[j]];
                adj[v[j]][v[last]] += adj[v[prev]][v[j]];
            }
            
            v.erase(v.begin() + prev);
            sz--;
        }
        
        long long min_reveals = n - 1;
        long long max_reveals = total_pairs - min_c + 1;
        cout << min_reveals << "\n" << max_reveals << "\n";
    }
    
    return 0;
}