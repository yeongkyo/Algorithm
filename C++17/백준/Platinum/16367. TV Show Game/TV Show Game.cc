#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

int K, N;
vector<int> adj[10005];
int dfn[10005], low[10005], scc[10005];
bool in_st[10005];
stack<int> st;
int timer, scc_cnt;

void tarjan(int u) {
    dfn[u] = low[u] = ++timer;
    st.push(u);
    in_st[u] = true;

    for (int v : adj[u]) {
        if (!dfn[v]) {
            tarjan(v);
            low[u] = min(low[u], low[v]);
        } else if (in_st[v]) {
            low[u] = min(low[u], dfn[v]);
        }
    }

    if (low[u] == dfn[u]) {
        scc_cnt++;
        while (true) {
            int t = st.top();
            st.pop();
            in_st[t] = false;
            scc[t] = scc_cnt;
            if (u == t) break;
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> K >> N)) return 0;

    for (int i = 0; i < N; ++i) {
        int l[3];
        char c[3];
        int lit[3], nlit[3];
        for (int j = 0; j < 3; ++j) {
            cin >> l[j] >> c[j];
            if (c[j] == 'R') {
                lit[j] = 2 * l[j];
                nlit[j] = 2 * l[j] + 1;
            } else {
                lit[j] = 2 * l[j] + 1;
                nlit[j] = 2 * l[j];
            }
        }

        for (int j = 0; j < 3; ++j) {
            int u = nlit[j];
            int v1 = lit[(j + 1) % 3];
            int v2 = lit[(j + 2) % 3];
            adj[u].push_back(v1);
            adj[nlit[(j + 1) % 3]].push_back(lit[j]);
            adj[u].push_back(v2);
            adj[nlit[(j + 2) % 3]].push_back(lit[j]);
        }
    }

    for (int i = 2; i <= 2 * K + 1; ++i) {
        if (!dfn[i]) tarjan(i);
    }

    for (int i = 1; i <= K; ++i) {
        if (scc[2 * i] == scc[2 * i + 1]) {
            cout << -1 << endl;
            return 0;
        }
    }

    string res = "";
    for (int i = 1; i <= K; ++i) {
        if (scc[2 * i] < scc[2 * i + 1]) res += 'R';
        else res += 'B';
    }
    cout << res << endl;

    return 0;
}