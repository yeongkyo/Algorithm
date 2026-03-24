#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const long long INF_COST = 1e18;
const int INF_FLOW = 1e9;

struct Edge {
    int to;
    int cap;
    int flow;
    long long cost;
    int rev;
};

vector<vector<Edge>> adj;

void add_edge(int from, int to, int cap, long long cost) {
    adj[from].push_back({to, cap, 0, cost, (int)adj[to].size()});
    adj[to].push_back({from, 0, 0, -cost, (int)adj[from].size() - 1});
}

void spfa(int S, int T, vector<long long>& pot) {
    int n = adj.size();
    pot.assign(n, INF_COST);
    pot[S] = 0;
    queue<int> q;
    vector<bool> in_q(n, false);
    q.push(S);
    in_q[S] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        in_q[u] = false;

        for (auto& e : adj[u]) {
            if (e.cap - e.flow > 0 && pot[e.to] > pot[u] + e.cost) {
                pot[e.to] = pot[u] + e.cost;
                if (!in_q[e.to]) {
                    q.push(e.to);
                    in_q[e.to] = true;
                }
            }
        }
    }
}

pair<int, long long> mcmf(int S, int T) {
    int n = adj.size();
    vector<long long> pot(n, 0);
    spfa(S, T, pot);

    int max_flow = 0;
    long long min_cost = 0;

    while (true) {
        vector<long long> dist(n, INF_COST);
        vector<int> parent_edge(n, -1);
        vector<int> parent_node(n, -1);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;

        dist[S] = 0;
        pq.push({0, S});

        while (!pq.empty()) {
            pair<long long, int> top_node = pq.top();
            pq.pop();
            long long d = top_node.first;
            int u = top_node.second;

            if (dist[u] < d) continue;

            for (int i = 0; i < adj[u].size(); i++) {
                Edge& e = adj[u][i];
                if (e.cap - e.flow > 0) {
                    long long nc = dist[u] + e.cost + pot[u] - pot[e.to];
                    if (dist[e.to] > nc) {
                        dist[e.to] = nc;
                        parent_node[e.to] = u;
                        parent_edge[e.to] = i;
                        pq.push({dist[e.to], e.to});
                    }
                }
            }
        }

        if (dist[T] == INF_COST) break;

        for (int i = 0; i < n; i++) {
            if (dist[i] != INF_COST) pot[i] += dist[i];
        }

        int push = INF_FLOW;
        int curr = T;
        while (curr != S) {
            int p = parent_node[curr];
            int idx = parent_edge[curr];
            push = min(push, adj[p][idx].cap - adj[p][idx].flow);
            curr = p;
        }

        curr = T;
        while (curr != S) {
            int p = parent_node[curr];
            int idx = parent_edge[curr];
            int rev_idx = adj[p][idx].rev;
            adj[p][idx].flow += push;
            adj[curr][rev_idx].flow -= push;
            min_cost += (long long)push * adj[p][idx].cost;
            curr = p;
        }

        max_flow += push;
    }

    return {max_flow, min_cost};
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<long long> A(N + 1), H(N + 1), L(N + 1);
    int root = 1;
    for (int i = 1; i <= N; ++i) {
        cin >> A[i];
        if (A[i] > A[root]) root = i;
    }
    for (int i = 1; i <= N; ++i) cin >> H[i];
    for (int i = 1; i <= N; ++i) cin >> L[i];

    int S = 0;
    int T = 2 * N + 1;
    adj.assign(T + 1, vector<Edge>());

    for (int i = 1; i <= N; ++i) {
        if (i != root) {
            add_edge(S, i, 1, 0);
            add_edge(i + N, T, L[i] - 1, 0);
        } else {
            add_edge(i + N, T, L[i], 0);
        }
    }

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            if (A[i] < A[j]) {
                long long weight = (A[i] ^ A[j]) - H[i] - H[j];
                add_edge(i, j + N, 1, -weight);
            }
        }
    }

    pair<int, long long> result = mcmf(S, T);
    cout << -result.second << "\n";

    return 0;
}