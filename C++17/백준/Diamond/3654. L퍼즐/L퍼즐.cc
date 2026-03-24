#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int n, m;
vector<string> grid;
vector<vector<int>> B_idx;

int v_count;
vector<vector<int>> adj;
vector<int> dfn, low, scc;
vector<bool> in_stack;
vector<int> st;
int timer, scc_cnt;

void add_edge(int u, int v) {
    adj[u].push_back(v);
}

void add_clause(int u, int v) {
    add_edge(u ^ 1, v);
    add_edge(v ^ 1, u);
}

void dfs(int u) {
    dfn[u] = low[u] = ++timer;
    st.push_back(u);
    in_stack[u] = true;

    for (int v : adj[u]) {
        if (!dfn[v]) {
            dfs(v);
            low[u] = min(low[u], low[v]);
        } else if (in_stack[v]) {
            low[u] = min(low[u], dfn[v]);
        }
    }

    if (low[u] == dfn[u]) {
        scc_cnt++;
        while (true) {
            int v = st.back();
            st.pop_back();
            in_stack[v] = false;
            scc[v] = scc_cnt;
            if (u == v) break;
        }
    }
}

void solve() {
    cin >> n >> m;
    grid.assign(n, "");
    B_idx.assign(n, vector<int>(m, -1));

    int b_cnt = 0;
    int w_cnt = 0;

    for (int i = 0; i < n; i++) {
        cin >> grid[i];
        for (int j = 0; j < m; j++) {
            if (grid[i][j] == 'B') {
                B_idx[i][j] = b_cnt++;
            } else if (grid[i][j] == 'W') {
                w_cnt++;
            }
        }
    }

    if (b_cnt * 2 != w_cnt) {
        cout << "NO\n";
        return;
    }

    if (b_cnt == 0) {
        cout << "YES\n";
        return;
    }

    int num_vars = 2 * b_cnt;
    adj.assign(2 * num_vars, vector<int>());

    auto X = [&](int i) { return i; };
    auto Y = [&](int i) { return i + b_cnt; };

    auto is_W = [&](int r, int c) {
        return r >= 0 && r < n && c >= 0 && c < m && grid[r][c] == 'W';
    };

    for (int r = 0; r < n; r++) {
        for (int c = 0; c < m; c++) {
            if (grid[r][c] == 'B') {
                int id = B_idx[r][c];

                if (!is_W(r, c - 1)) add_clause(2 * X(id) + 1, 2 * X(id) + 1);
                if (!is_W(r, c + 1)) add_clause(2 * X(id), 2 * X(id));
                if (!is_W(r - 1, c)) add_clause(2 * Y(id) + 1, 2 * Y(id) + 1);
                if (!is_W(r + 1, c)) add_clause(2 * Y(id), 2 * Y(id));
            } else if (grid[r][c] == 'W') {
                vector<int> lits;
                if (c - 1 >= 0 && grid[r][c - 1] == 'B') lits.push_back(2 * X(B_idx[r][c - 1]) + 1);
                if (c + 1 < m && grid[r][c + 1] == 'B') lits.push_back(2 * X(B_idx[r][c + 1]));
                if (r - 1 >= 0 && grid[r - 1][c] == 'B') lits.push_back(2 * Y(B_idx[r - 1][c]) + 1);
                if (r + 1 < n && grid[r + 1][c] == 'B') lits.push_back(2 * Y(B_idx[r + 1][c]));

                for (size_t i = 0; i < lits.size(); i++) {
                    for (size_t j = i + 1; j < lits.size(); j++) {
                        add_clause(lits[i] ^ 1, lits[j] ^ 1);
                    }
                }
            }
        }
    }

    dfn.assign(2 * num_vars, 0);
    low.assign(2 * num_vars, 0);
    scc.assign(2 * num_vars, 0);
    in_stack.assign(2 * num_vars, false);
    timer = scc_cnt = 0;
    st.clear();

    for (int i = 0; i < 2 * num_vars; i++) {
        if (!dfn[i]) dfs(i);
    }

    bool possible = true;
    for (int i = 0; i < num_vars; i++) {
        if (scc[2 * i] == scc[2 * i + 1]) {
            possible = false;
            break;
        }
    }

    if (possible) cout << "YES\n";
    else cout << "NO\n";
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