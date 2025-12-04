#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const int INF = 1e9;
const int MAX_V = 205;

struct Edge {
    int to, cap, flow, cost, rev;
};

vector<Edge> adj[MAX_V];
int dist[MAX_V], p[MAX_V], p_edge[MAX_V];
bool in_queue[MAX_V];

void add_edge(int u, int v, int cap, int cost) {
    adj[u].push_back({v, cap, 0, cost, (int)adj[v].size()});
    adj[v].push_back({u, 0, 0, -cost, (int)adj[u].size() - 1});
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M;
    cin >> N >> M;

    int S = 0;
    int T = N + M + 1;

    vector<int> A(N + 1);
    for (int i = 1; i <= N; i++) {
        cin >> A[i];
        add_edge(M + i, T, A[i], 0);
    }

    vector<int> B(M + 1);
    for (int i = 1; i <= M; i++) {
        cin >> B[i];
        add_edge(S, i, B[i], 0);
    }

    vector<vector<int>> C(M + 1, vector<int>(N + 1));
    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            cin >> C[i][j];
        }
    }

    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            int d;
            cin >> d;
            if (C[i][j] > 0) {
                add_edge(i, M + j, C[i][j], d);
            }
        }
    }

    int total_flow = 0;
    int total_cost = 0;

    while (true) {
        fill(dist, dist + MAX_V, INF);
        fill(in_queue, in_queue + MAX_V, false);
        fill(p, p + MAX_V, -1);
        fill(p_edge, p_edge + MAX_V, -1);

        queue<int> q;
        dist[S] = 0;
        in_queue[S] = true;
        q.push(S);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            in_queue[u] = false;

            for (int i = 0; i < adj[u].size(); i++) {
                Edge &e = adj[u][i];
                if (e.cap - e.flow > 0 && dist[e.to] > dist[u] + e.cost) {
                    dist[e.to] = dist[u] + e.cost;
                    p[e.to] = u;
                    p_edge[e.to] = i;
                    if (!in_queue[e.to]) {
                        q.push(e.to);
                        in_queue[e.to] = true;
                    }
                }
            }
        }

        if (dist[T] == INF) break;

        int flow = INF;
        for (int i = T; i != S; i = p[i]) {
            int prev = p[i];
            int idx = p_edge[i];
            flow = min(flow, adj[prev][idx].cap - adj[prev][idx].flow);
        }

        for (int i = T; i != S; i = p[i]) {
            int prev = p[i];
            int idx = p_edge[i];
            adj[prev][idx].flow += flow;
            int rev_idx = adj[prev][idx].rev;
            adj[i][rev_idx].flow -= flow;
            total_cost += flow * adj[prev][idx].cost;
        }

        total_flow += flow;
    }

    cout << total_flow << "\n";
    cout << total_cost << "\n";

    return 0;
}