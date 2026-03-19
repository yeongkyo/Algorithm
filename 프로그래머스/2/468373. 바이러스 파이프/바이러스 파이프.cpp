#include <vector>
#include <bitset>
#include <queue>
#include <algorithm>

using namespace std;

int max_inf = 0;
int N;
bitset<105> reach[3][105];

void dfs(int k_left, int last_type, bitset<105> S) {
    max_inf = max(max_inf, (int)S.count());
    if (k_left == 0) return;
    
    for (int t = 0; t < 3; ++t) {
        if (t == last_type) continue;
        
        bitset<105> next_S = S;
        for (int i = 1; i <= N; ++i) {
            if (S[i]) {
                next_S |= reach[t][i];
            }
        }
        
        if (next_S != S) {
            dfs(k_left - 1, t, next_S);
        }
    }
}

int solution(int n, int infection, vector<vector<int>> edges, int k) {
    N = n;
    max_inf = 0;
    
    vector<int> adj[3][105];
    for (auto& e : edges) {
        int u = e[0];
        int v = e[1];
        int type = e[2] - 1;
        adj[type][u].push_back(v);
        adj[type][v].push_back(u);
    }
    
    for (int t = 0; t < 3; ++t) {
        for (int i = 1; i <= n; ++i) {
            reach[t][i].reset();
            queue<int> q;
            q.push(i);
            vector<bool> vis(n + 1, false);
            vis[i] = true;
            reach[t][i].set(i);
            
            while (!q.empty()) {
                int curr = q.front();
                q.pop();
                for (int nxt : adj[t][curr]) {
                    if (!vis[nxt]) {
                        vis[nxt] = true;
                        reach[t][i].set(nxt);
                        q.push(nxt);
                    }
                }
            }
        }
    }
    
    bitset<105> start_S;
    start_S.set(infection);
    
    dfs(k, -1, start_S);
    
    return max_inf;
}