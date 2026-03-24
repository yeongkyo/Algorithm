#pragma GCC optimize("O3")
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <set>
#include <string>

using namespace std;

const __int128 INF = (__int128)1 << 100;

struct Edge {
    int to;
    int rev;
    __int128 cap;
    __int128 flow;
};

struct CapUpdate {
    int u;
    int idx;
    long long C;
    long long R;
};

int N, M, P;
long long K;
vector<vector<Edge>> adj;
vector<CapUpdate> updates;
vector<int> level_arr;
vector<int> ptr;

int S, T;

long long C_cap[305][305];
long long R_cap[305][305];

struct Line {
    __int128 A;
    __int128 B;
};

vector<Line> lines;

void add_line(__int128 A, __int128 B) {
    for (auto& l : lines) {
        if (l.A == A && l.B == B) return;
    }
    lines.push_back({A, B});
}

set<long long> evaluated;
__int128 max_F_found = -1;

void add_edge(int from, int to, __int128 cap, bool is_internal = false, long long C = 0, long long R = 0) {
    int idx_from = adj[from].size();
    int idx_to = adj[to].size();
    adj[from].push_back({to, idx_to, cap, 0});
    adj[to].push_back({from, idx_from, 0, 0});
    if (is_internal) {
        updates.push_back({from, idx_from, C, R});
    }
}

bool bfs() {
    fill(level_arr.begin(), level_arr.end(), -1);
    level_arr[S] = 0;
    queue<int> q;
    q.push(S);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (auto& edge : adj[v]) {
            if (edge.cap - edge.flow > 0 && level_arr[edge.to] == -1) {
                level_arr[edge.to] = level_arr[v] + 1;
                q.push(edge.to);
            }
        }
    }
    return level_arr[T] != -1;
}

__int128 dfs(int v, __int128 pushed) {
    if (pushed == 0) return 0;
    if (v == T) return pushed;
    for (int& cid = ptr[v]; cid < adj[v].size(); ++cid) {
        auto& edge = adj[v][cid];
        int tr = edge.to;
        if (level_arr[v] + 1 != level_arr[tr] || edge.cap - edge.flow == 0) continue;
        __int128 push = dfs(tr, min(pushed, edge.cap - edge.flow));
        if (push == 0) continue;
        edge.flow += push;
        adj[tr][edge.rev].flow -= push;
        return push;
    }
    return 0;
}

__int128 get_flow(long long t) {
    for (auto& up : updates) {
        adj[up.u][up.idx].cap = up.C + (__int128)t * up.R;
    }
    for (int i = 0; i < adj.size(); ++i) {
        for (auto& e : adj[i]) {
            e.flow = 0;
        }
    }
    __int128 flow = 0;
    while (bfs()) {
        fill(ptr.begin(), ptr.end(), 0);
        while (true) {
            __int128 pushed = dfs(S, INF);
            if (!pushed) break;
            flow += pushed;
        }
    }
    return flow;
}

int get_id(int r, int c) {
    return (r - 1) * M + c;
}

void eval(long long t) {
    if (evaluated.count(t)) return;
    __int128 flow = get_flow(t);
    __int128 A = 0, B = 0;
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            int in_node = get_id(i, j);
            int out_node = in_node + N * M;
            if (level_arr[in_node] != -1 && level_arr[out_node] == -1) {
                A += R_cap[i][j];
                B += C_cap[i][j];
            }
        }
    }
    add_line(A, B);
    if (flow > max_F_found) max_F_found = flow;
    evaluated.insert(t);
}

void print(__int128 x) {
    if (x == 0) {
        cout << 0 << "\n";
        return;
    }
    string s;
    while (x > 0) {
        s += (char)('0' + (int)(x % 10));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    cout << s << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    if (!(cin >> N >> M >> P >> K)) return 0;

    int total_nodes = 2 * N * M + 2;
    adj.resize(total_nodes);
    level_arr.resize(total_nodes);
    ptr.resize(total_nodes);

    S = 0;
    T = total_nodes - 1;

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            cin >> C_cap[i][j];
        }
    }

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            cin >> R_cap[i][j];
        }
    }

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            int in_node = get_id(i, j);
            int out_node = in_node + N * M;
            add_edge(in_node, out_node, 0, true, C_cap[i][j], R_cap[i][j]);
        }
    }

    add_edge(S, get_id(1, 1), INF);
    add_edge(get_id(N, M) + N * M, T, INF);

    for (int i = 0; i < P; ++i) {
        int x1, y1, x2, y2;
        cin >> x1 >> y1 >> x2 >> y2;
        int u_out = get_id(x1, y1) + N * M;
        int v_in = get_id(x2, y2);
        add_edge(u_out, v_in, INF);

        int v_out = get_id(x2, y2) + N * M;
        int u_in = get_id(x1, y1);
        add_edge(v_out, u_in, INF);
    }

    eval(0);
    if (K != 0) eval(K);

    while (true) {
        vector<long long> candidates;
        candidates.push_back(0);
        candidates.push_back(K);
        
        for (size_t i = 0; i < lines.size(); ++i) {
            for (size_t j = i + 1; j < lines.size(); ++j) {
                __int128 A1 = lines[i].A;
                __int128 B1 = lines[i].B;
                __int128 A2 = lines[j].A;
                __int128 B2 = lines[j].B;
                if (A1 == A2) continue;
                
                __int128 num = B2 - B1;
                __int128 den = A1 - A2;
                if (den < 0) {
                    num = -num;
                    den = -den;
                }
                
                __int128 t_floor = num / den;
                if (num % den < 0) t_floor--;
                
                __int128 t_ceil = t_floor + (num % den != 0 ? 1 : 0);
                
                if (t_floor >= 0 && t_floor <= K) candidates.push_back((long long)t_floor);
                if (t_ceil >= 0 && t_ceil <= K) candidates.push_back((long long)t_ceil);
            }
        }
        
        __int128 best_U = -1;
        long long best_t = -1;
        
        for (long long t : candidates) {
            __int128 U_t = -1;
            bool first = true;
            for (auto& l : lines) {
                __int128 val = l.A * t + l.B;
                if (first || val < U_t) {
                    U_t = val;
                    first = false;
                }
            }
            
            if (U_t > best_U) {
                best_U = U_t;
                best_t = t;
            } else if (U_t == best_U) {
                if (evaluated.find(t) == evaluated.end() && evaluated.find(best_t) != evaluated.end()) {
                    best_t = t;
                }
            }
        }
        
        if (best_U <= max_F_found) {
            break;
        }
        
        if (evaluated.find(best_t) != evaluated.end()) {
            break; 
        }
        
        eval(best_t);
    }

    print(max_F_found);

    return 0;
}