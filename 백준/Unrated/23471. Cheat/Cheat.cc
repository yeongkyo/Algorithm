#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void solve() {
    int n, m;
    if (!(cin >> n >> m)) return;

    int L = 1, R = n;
    vector<int> max_u_at(n + 1, 0);
    vector<pair<int, int>> forward_edges;

    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        if (u >= v) {
            L = max(L, v);
            R = min(R, u);
            max_u_at[v] = max(max_u_at[v], u);
        } else {
            forward_edges.push_back({u, v});
        }
    }

    if (L > R) {
        cout << 0 << "\n";
        return;
    }

    vector<int> max_u(n + 1, 0);
    int current_max = 0;
    for (int i = 1; i <= n; ++i) {
        current_max = max(current_max, max_u_at[i]);
        max_u[i] = current_max;
    }

    vector<int> diff(n + 2, 0);
    for (auto edge : forward_edges) {
        int a = edge.first;
        int b = edge.second;
        
        if (max_u[a] >= b) {
            if (a + 1 <= b - 1) {
                diff[a + 1]++;
                diff[b]--;
            }
        }
    }

    vector<int> ans;
    int current_invalid = 0;
    
    for (int i = 1; i <= n; ++i) {
        current_invalid += diff[i];
        if (i >= L && i <= R && current_invalid == 0) {
            ans.push_back(i);
        }
    }

    cout << ans.size();
    for (int x : ans) {
        cout << " " << x;
    }
    cout << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int z;
    if (cin >> z) {
        while (z--) {
            solve();
        }
    }
    return 0;
}