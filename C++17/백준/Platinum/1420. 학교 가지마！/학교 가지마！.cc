#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

const int INF = 1e9;
int N, M;
string grid[101];
int dr[] = { -1, 1, 0, 0 };
int dc[] = { 0, 0, -1, 1 };

struct Edge {
    int to, capacity, flow, rev;
};

vector<Edge> adj[20005];
int level[20005];
int iter[20005];

void add_edge(int from, int to, int cap) {
    adj[from].push_back({ to, cap, 0, (int)adj[to].size() });
    adj[to].push_back({ from, 0, 0, (int)adj[from].size() - 1 });
}

bool bfs(int s, int t) {
    fill(level, level + 20005, -1);
    level[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (auto& e : adj[v]) {
            if (e.capacity - e.flow > 0 && level[e.to] < 0) {
                level[e.to] = level[v] + 1;
                q.push(e.to);
            }
        }
    }
    return level[t] != -1;
}

int dfs(int v, int t, int f) {
    if (v == t) return f;
    for (int& i = iter[v]; i < adj[v].size(); ++i) {
        Edge& e = adj[v][i];
        if (e.capacity - e.flow > 0 && level[v] < level[e.to]) {
            int d = dfs(e.to, t, min(f, e.capacity - e.flow));
            if (d > 0) {
                e.flow += d;
                adj[e.to][e.rev].flow -= d;
                return d;
            }
        }
    }
    return 0;
}

int max_flow(int s, int t) {
    int flow = 0;
    while (bfs(s, t)) {
        fill(iter, iter + 20005, 0);
        int f;
        while ((f = dfs(s, t, INF)) > 0) {
            flow += f;
        }
    }
    return flow;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> N >> M;
    int sr, sc, tr, tc;
    for (int i = 0; i < N; ++i) {
        cin >> grid[i];
        for (int j = 0; j < M; ++j) {
            if (grid[i][j] == 'K') { sr = i; sc = j; }
            if (grid[i][j] == 'H') { tr = i; tc = j; }
        }
    }

    for (int i = 0; i < 4; ++i) {
        int nr = sr + dr[i];
        int nc = sc + dc[i];
        if (nr >= 0 && nr < N && nc >= 0 && nc < M && grid[nr][nc] == 'H') {
            cout << -1 << endl;
            return 0;
        }
    }

    int source = (sr * M + sc) * 2 + 1;
    int sink = (tr * M + tc) * 2;

    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            if (grid[r][c] == '#') continue;

            int in = (r * M + c) * 2;
            int out = (r * M + c) * 2 + 1;
            int cap = (grid[r][c] == 'K' || grid[r][c] == 'H') ? INF : 1;
            add_edge(in, out, cap);

            for (int i = 0; i < 4; ++i) {
                int nr = r + dr[i];
                int nc = c + dc[i];
                if (nr >= 0 && nr < N && nc >= 0 && nc < M && grid[nr][nc] != '#') {
                    int n_in = (nr * M + nc) * 2;
                    add_edge(out, n_in, INF);
                }
            }
        }
    }

    int result = max_flow(source, sink);
    cout << (result >= INF ? -1 : result) << endl;

    return 0;
}