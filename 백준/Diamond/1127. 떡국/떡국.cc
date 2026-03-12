#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

const int INF = 1e9;

struct Edge {
    int to;
    int cap;
    int flow;
    int rev;
};

vector<vector<Edge>> adj_flow;
vector<int> level;
vector<int> ptr;

void add_directed_edge(int from, int to, int cap) {
    adj_flow[from].push_back({to, cap, 0, (int)adj_flow[to].size()});
    adj_flow[to].push_back({from, 0, 0, (int)adj_flow[from].size() - 1});
}

void add_undirected_edge(int u, int v, int cap) {
    adj_flow[u].push_back({v, cap, 0, (int)adj_flow[v].size()});
    adj_flow[v].push_back({u, cap, 0, (int)adj_flow[u].size() - 1});
}

bool bfs(int s, int t) {
    fill(level.begin(), level.end(), -1);
    level[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (auto& edge : adj_flow[v]) {
            if (edge.cap - edge.flow > 0 && level[edge.to] == -1) {
                level[edge.to] = level[v] + 1;
                q.push(edge.to);
            }
        }
    }
    return level[t] != -1;
}

int dfs(int v, int t, int pushed) {
    if (pushed == 0) return 0;
    if (v == t) return pushed;
    for (int& cid = ptr[v]; cid < adj_flow[v].size(); ++cid) {
        auto& edge = adj_flow[v][cid];
        int tr = edge.to;
        if (level[v] + 1 != level[tr] || edge.cap - edge.flow == 0) continue;
        int push = dfs(tr, t, min(pushed, edge.cap - edge.flow));
        if (push == 0) continue;
        edge.flow += push;
        adj_flow[tr][edge.rev].flow -= push;
        return push;
    }
    return 0;
}

int dinic(int s, int t) {
    int flow = 0;
    while (bfs(s, t)) {
        fill(ptr.begin(), ptr.end(), 0);
        while (int pushed = dfs(s, t, INF)) {
            flow += pushed;
        }
    }
    return flow;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<string> adj(N);
    for (int i = 0; i < N; ++i) {
        cin >> adj[i];
    }

    int K;
    cin >> K;

    vector<vector<bool>> has_guard(N, vector<bool>(K, false));
    for (int i = 0; i < N; ++i) {
        int cnt;
        cin >> cnt;
        for (int j = 0; j < cnt; ++j) {
            int c;
            cin >> c;
            has_guard[i][c] = true;
        }
    }

    vector<vector<bool>> has_office(N, vector<bool>(K, false));
    for (int i = 0; i < N; ++i) {
        int cnt;
        cin >> cnt;
        for (int j = 0; j < cnt; ++j) {
            int c;
            cin >> c;
            has_office[i][c] = true;
        }
    }

    int total_conflicts = 0;
    int S = N, T = N + 1;

    for (int c = 0; c < K; ++c) {
        adj_flow.assign(N + 2, vector<Edge>());
        level.resize(N + 2);
        ptr.resize(N + 2);

        for (int i = 0; i < N; ++i) {
            if (has_guard[i][c]) {
                add_directed_edge(S, i, INF);
            } else if (!has_office[i][c]) {
                add_directed_edge(i, T, INF);
            }
        }

        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (adj[i][j] == '1') {
                    add_undirected_edge(i, j, 1);
                }
            }
        }

        total_conflicts += dinic(S, T);
    }

    cout << total_conflicts << "\n";

    return 0;
}