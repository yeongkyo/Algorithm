#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <stack>
#include <algorithm>

using namespace std;

enum InstType { NCS, CS, SET, TEST };

struct Inst {
    InstType type;
    int var;
    int val;
    int next0;
    int next1;
};

struct State {
    int ip0, ip1, A, B, C;
    int id() const {
        return ip0 * 80 + ip1 * 8 + A * 4 + B * 2 + C;
    }
    static State decode(int id) {
        State s;
        s.C = id % 2; id /= 2;
        s.B = id % 2; id /= 2;
        s.A = id % 2; id /= 2;
        s.ip1 = id % 10; id /= 10;
        s.ip0 = id;
        return s;
    }
};

bool check_violation(int mode, const vector<Inst>& code0, const vector<Inst>& code1, const bool reach[800], const int adj[800][2]) {
    vector<vector<pair<int, int>>> g(800);
    for (int u = 0; u < 800; u++) {
        if (!reach[u]) continue;
        State s = State::decode(u);
        
        bool add0 = true;
        if (mode == 0 && code0[s.ip0].type == CS) add0 = false;
        if (mode == 1 && code0[s.ip0].type == CS) add0 = false;
        if (add0) g[u].push_back({adj[u][0], 0});

        bool add1 = true;
        if (mode == 0 && code1[s.ip1].type == CS) add1 = false;
        if (mode == 2 && code1[s.ip1].type == CS) add1 = false;
        if (add1) g[u].push_back({adj[u][1], 1});
    }

    int dfn[800] = {0}, low[800] = {0}, timer = 0;
    bool in_stack[800] = {false};
    stack<int> st;
    int scc_id[800] = {0};
    int scc_count = 0;
    vector<vector<int>> sccs;

    auto dfs = [&](auto& self, int u) -> void {
        dfn[u] = low[u] = ++timer;
        st.push(u);
        in_stack[u] = true;
        for (auto& edge : g[u]) {
            int v = edge.first;
            if (!dfn[v]) {
                self(self, v);
                low[u] = min(low[u], low[v]);
            } else if (in_stack[v]) {
                low[u] = min(low[u], dfn[v]);
            }
        }
        if (low[u] == dfn[u]) {
            scc_count++;
            vector<int> comp;
            while (true) {
                int v = st.top(); st.pop();
                in_stack[v] = false;
                scc_id[v] = scc_count;
                comp.push_back(v);
                if (u == v) break;
            }
            sccs.push_back(comp);
        }
    };

    for (int i = 0; i < 800; i++) {
        if (reach[i] && !dfn[i]) dfs(dfs, i);
    }

    for (const auto& comp : sccs) {
        bool has_T0 = false, has_T1 = false;
        bool all_ip0_NCS = true, all_ip1_NCS = true;
        bool has_edges = false;

        for (int u : comp) {
            State s = State::decode(u);
            if (code0[s.ip0].type != NCS) all_ip0_NCS = false;
            if (code1[s.ip1].type != NCS) all_ip1_NCS = false;

            for (auto& edge : g[u]) {
                int v = edge.first;
                int thread = edge.second;
                if (scc_id[v] == scc_id[u]) {
                    has_edges = true;
                    if (thread == 0) has_T0 = true;
                    if (thread == 1) has_T1 = true;
                }
            }
        }

        if (!has_edges) continue;

        if (mode == 0) { // Deadlock
            if ((has_T0 && has_T1) || (has_T0 && !has_T1 && all_ip1_NCS) || (!has_T0 && has_T1 && all_ip0_NCS)) {
                return true; 
            }
        } else if (mode == 1) { // Starvation T0
            if (has_T0 && (has_T1 || (!has_T1 && all_ip1_NCS))) {
                return true; 
            }
        } else if (mode == 2) { // Starvation T1
            if (has_T1 && (has_T0 || (!has_T0 && all_ip0_NCS))) {
                return true;
            }
        }
    }
    return false;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int m1, m2;
    while (cin >> m1 >> m2) {
        vector<Inst> code0(m1 + 1), code1(m2 + 1);

        auto read_code = [](int m, vector<Inst>& code) {
            for (int i = 0; i < m; i++) {
                int lino;
                string op;
                cin >> lino >> op;
                Inst inst;
                if (op == "NCS" || op == "CS") {
                    inst.type = (op == "NCS") ? NCS : CS;
                    cin >> inst.next0;
                } else if (op == "SET") {
                    inst.type = SET;
                    string var;
                    cin >> var >> inst.val >> inst.next0;
                    inst.var = var[0] - 'A';
                } else if (op == "TEST") {
                    inst.type = TEST;
                    string var;
                    cin >> var >> inst.next0 >> inst.next1;
                    inst.var = var[0] - 'A';
                }
                code[lino] = inst;
            }
        };

        read_code(m1, code0);
        read_code(m2, code1);

        int adj[800][2] = {0};
        for (int i = 0; i < 800; i++) {
            State s = State::decode(i);
            if (s.ip0 < 1 || s.ip0 > m1 || s.ip1 < 1 || s.ip1 > m2) continue;

            auto get_next = [&](const Inst& inst, int& n_ip, int& n_A, int& n_B, int& n_C) {
                if (inst.type == NCS || inst.type == CS) {
                    n_ip = inst.next0;
                } else if (inst.type == SET) {
                    if (inst.var == 0) n_A = inst.val;
                    else if (inst.var == 1) n_B = inst.val;
                    else n_C = inst.val;
                    n_ip = inst.next0;
                } else if (inst.type == TEST) {
                    int v = (inst.var == 0) ? s.A : (inst.var == 1 ? s.B : s.C);
                    n_ip = (v == 0) ? inst.next0 : inst.next1;
                }
            };

            int n_ip0 = s.ip0, n_A0 = s.A, n_B0 = s.B, n_C0 = s.C;
            get_next(code0[s.ip0], n_ip0, n_A0, n_B0, n_C0);
            adj[i][0] = State{n_ip0, s.ip1, n_A0, n_B0, n_C0}.id();

            int n_ip1 = s.ip1, n_A1 = s.A, n_B1 = s.B, n_C1 = s.C;
            get_next(code1[s.ip1], n_ip1, n_A1, n_B1, n_C1);
            adj[i][1] = State{s.ip0, n_ip1, n_A1, n_B1, n_C1}.id();
        }

        bool reach[800] = {false};
        queue<int> q;
        int start_id = State{1, 1, 0, 0, 0}.id();
        q.push(start_id);
        reach[start_id] = true;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (!reach[adj[u][0]]) { reach[adj[u][0]] = true; q.push(adj[u][0]); }
            if (!reach[adj[u][1]]) { reach[adj[u][1]] = true; q.push(adj[u][1]); }
        }

        bool me_satisfied = true;
        for (int i = 0; i < 800; i++) {
            if (reach[i]) {
                State s = State::decode(i);
                if (code0[s.ip0].type == CS && code1[s.ip1].type == CS) {
                    me_satisfied = false;
                    break;
                }
            }
        }

        bool df_violated = check_violation(0, code0, code1, reach, adj);
        bool sf0_violated = check_violation(1, code0, code1, reach, adj);
        bool sf1_violated = check_violation(2, code0, code1, reach, adj);

        string ans = "";
        ans += me_satisfied ? "Y" : "N";
        ans += !df_violated ? "Y" : "N";
        ans += (!sf0_violated && !sf1_violated) ? "Y" : "N";
        
        cout << ans << "\n";
    }

    return 0;
}