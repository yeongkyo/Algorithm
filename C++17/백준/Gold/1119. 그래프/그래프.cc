#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<string> adj(n);
    for (int i = 0; i < n; ++i) {
        cin >> adj[i];
    }

    if (n == 1) {
        cout << 0 << "\n";
        return 0;
    }

    int edges2 = 0;
    for (int i = 0; i < n; ++i) {
        int deg = 0;
        for (int j = 0; j < n; ++j) {
            if (adj[i][j] == 'Y') {
                deg++;
                edges2++;
            }
        }
        if (deg == 0) {
            cout << -1 << "\n";
            return 0;
        }
    }

    int m = edges2 / 2;
    if (m < n - 1) {
        cout << -1 << "\n";
        return 0;
    }

    vector<bool> visited(n, false);
    int comps = 0;

    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            comps++;
            vector<int> st;
            st.push_back(i);
            visited[i] = true;
            while (!st.empty()) {
                int u = st.back();
                st.pop_back();
                for (int v = 0; v < n; ++v) {
                    if (adj[u][v] == 'Y' && !visited[v]) {
                        visited[v] = true;
                        st.push_back(v);
                    }
                }
            }
        }
    }

    cout << comps - 1 << "\n";

    return 0;
}