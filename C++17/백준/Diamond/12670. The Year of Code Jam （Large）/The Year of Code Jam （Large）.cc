#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

const long long INF = 1e18;

struct Edge {
    int to;
    long long cap, flow;
    int rev;
};

int N, M;
vector<vector<Edge>> adj;
vector<int> level;
vector<int> ptr;

void add_edge(int from, int to, long long cap) {
    adj[from].push_back({to, cap, 0, (int)adj[to].size()});
    adj[to].push_back({from, 0, 0, (int)adj[from].size() - 1});
}

bool bfs(int s, int t) {
    fill(level.begin(), level.end(), -1);
    level[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (auto& edge : adj[v]) {
            if (edge.cap - edge.flow > 0 && level[edge.to] == -1) {
                level[edge.to] = level[v] + 1;
                q.push(edge.to);
            }
        }
    }
    return level[t] != -1;
}

long long dfs(int v, int t, long long pushed) {
    if (pushed == 0) return 0;
    if (v == t) return pushed;
    for (int& cid = ptr[v]; cid < adj[v].size(); ++cid) {
        auto& edge = adj[v][cid];
        int tr = edge.to;
        if (level[v] + 1 != level[tr] || edge.cap - edge.flow == 0) continue;
        long long push = dfs(tr, t, min(pushed, edge.cap - edge.flow));
        if (push == 0) continue;
        edge.flow += push;
        adj[tr][edge.rev].flow -= push;
        return push;
    }
    return 0;
}

long long dinic(int s, int t) {
    long long flow = 0;
    while (bfs(s, t)) {
        fill(ptr.begin(), ptr.end(), 0);
        while (long long pushed = dfs(s, t, INF)) {
            flow += pushed;
        }
    }
    return flow;
}

void solve(int tc) {
    cin >> N >> M;
    vector<string> grid(N);
    int base_score = 0;
    for (int i = 0; i < N; ++i) {
        cin >> grid[i];
        for (int j = 0; j < M; ++j) {
            if (grid[i][j] == '#' || grid[i][j] == '?') {
                base_score += 4;
            }
        }
    }

    int S = N * M;
    int T = N * M + 1;
    adj.assign(T + 1, vector<Edge>());
    level.resize(T + 1);
    ptr.resize(T + 1);

    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};

    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            int u = r * M + c;
            char ch = grid[r][c];

            if ((r + c) % 2 == 0) {
                if (ch == '#') add_edge(S, u, INF);
                else if (ch == '?') add_edge(S, u, 4);
                else if (ch == '.') add_edge(u, T, INF);

                for (int i = 0; i < 4; ++i) {
                    int nr = r + dx[i];
                    int nc = c + dy[i];
                    if (nr >= 0 && nr < N && nc >= 0 && nc < M) {
                        int v = nr * M + nc;
                        add_edge(u, v, 2);
                    }
                }
            } else {
                if (ch == '#') add_edge(u, T, INF);
                else if (ch == '?') add_edge(u, T, 4);
                else if (ch == '.') add_edge(S, u, INF);
            }
        }
    }

    long long min_cut = dinic(S, T);
    cout << "Case #" << tc << ": " << base_score - min_cut << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int T;
    if (cin >> T) {
        for (int i = 1; i <= T; ++i) {
            solve(i);
        }
    }
    return 0;
}