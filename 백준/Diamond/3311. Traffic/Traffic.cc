#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>

using namespace std;

const int MAXN = 300005;
const int INF = 1e9 + 7;

struct Point {
    int id;
    int x, y;
};

int n, m;
int A_width, B_height;
Point points[MAXN];
vector<int> adj[MAXN];

// Tarjan's SCC variables
int dfn[MAXN], low[MAXN], scc_id[MAXN], scc_cnt = 0, timer = 0;
bool in_stack[MAXN];
stack<int> st;

// DAG variables
vector<int> dag_adj[MAXN];
bool scc_has_west[MAXN]; 
bool scc_reachable_from_west[MAXN]; // 서쪽에서 도달 가능한 SCC인지 여부
int scc_min[MAXN], scc_max[MAXN];

// Sorting comparators
bool compareYDesc(const Point& a, const Point& b) {
    if (a.y != b.y) return a.y > b.y;
    return a.x < b.x; 
}

void tarjan(int u) {
    dfn[u] = low[u] = ++timer;
    st.push(u);
    in_stack[u] = true;

    for (int v : adj[u]) {
        if (!dfn[v]) {
            tarjan(v);
            low[u] = min(low[u], low[v]);
        } else if (in_stack[v]) {
            low[u] = min(low[u], dfn[v]);
        }
    }

    if (low[u] == dfn[u]) {
        while (true) {
            int v = st.top();
            st.pop();
            in_stack[v] = false;
            scc_id[v] = scc_cnt;
            if (u == v) break;
        }
        scc_cnt++;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> n >> m >> A_width >> B_height;

    vector<Point> west_points;
    vector<Point> east_points;

    for (int i = 1; i <= n; i++) {
        points[i].id = i;
        cin >> points[i].x >> points[i].y;
        if (points[i].x == 0) west_points.push_back(points[i]);
        if (points[i].x == A_width) east_points.push_back(points[i]);
    }

    for (int i = 0; i < m; i++) {
        int u, v, k;
        cin >> u >> v >> k;
        adj[u].push_back(v);
        if (k == 2) adj[v].push_back(u);
    }

    // 1. SCC Decomposition
    for (int i = 1; i <= n; i++) {
        if (!dfn[i]) tarjan(i);
    }

    // 2. Build DAG
    for (int u = 1; u <= n; u++) {
        if (points[u].x == 0) scc_has_west[scc_id[u]] = true;
        
        for (int v : adj[u]) {
            if (scc_id[u] != scc_id[v]) {
                dag_adj[scc_id[u]].push_back(scc_id[v]);
            }
        }
    }

    // Remove duplicate edges in DAG
    for (int i = 0; i < scc_cnt; i++) {
        sort(dag_adj[i].begin(), dag_adj[i].end());
        dag_adj[i].erase(unique(dag_adj[i].begin(), dag_adj[i].end()), dag_adj[i].end());
    }

    // 3. Mark SCCs reachable from West (Top-Down on DAG)
    // Tarjan generates SCCs in reverse topological order (Sink -> Source).
    // So iterating from scc_cnt-1 down to 0 gives topological order (Source -> Sink).
    for (int i = scc_cnt - 1; i >= 0; i--) {
        if (scc_has_west[i]) scc_reachable_from_west[i] = true;
        if (scc_reachable_from_west[i]) {
            for (int next : dag_adj[i]) {
                scc_reachable_from_west[next] = true;
            }
        }
    }

    // 4. Identify reachable East nodes and assign ranks
    vector<Point> valid_east;
    for (const auto& p : east_points) {
        if (scc_reachable_from_west[scc_id[p.id]]) {
            valid_east.push_back(p);
        }
    }
    // Sort by Y descending
    sort(valid_east.begin(), valid_east.end(), compareYDesc);

    // Map East Node ID -> Rank (1 to K)
    // We only need to store min/max rank in each SCC
    for (int i = 0; i < scc_cnt; i++) {
        scc_min[i] = INF;
        scc_max[i] = -INF;
    }

    for (int i = 0; i < valid_east.size(); i++) {
        int sid = scc_id[valid_east[i].id];
        int rank = i + 1;
        scc_min[sid] = min(scc_min[sid], rank);
        scc_max[sid] = max(scc_max[sid], rank);
    }

    // 5. DP to propagate min/max ranges (Bottom-Up on DAG)
    // Iterate from 0 to scc_cnt-1 (Sink -> Source)
    for (int i = 0; i < scc_cnt; i++) {
        // We need to pull information from children? 
        // No, in Sink->Source order (0 to scc_cnt-1), when we visit U, 
        // its children V (which appear earlier in topological order, i.e., later in Tarjan loop?)
        // Wait. Tarjan logic:
        // dfs(u) calls dfs(v). v finishes first -> v added to SCC list first (lower ID).
        // So edges go from High ID to Low ID?
        // Let's check: U -> V. dfs(U) calls dfs(V).
        // If no cycle, V finishes, gets ID 0. Then U finishes, gets ID 1.
        // So Edge is 1 -> 0 (High to Low).
        // So to propagate from Children(0) to Parent(1), we iterate 0 to scc_cnt-1.
        // When we are at U(1), we look at neighbors. Neighbor V(0) is already processed.
        // Wait, normally DP uses value of *next* node.
        // If U -> V, min[u] = min(min[u], min[v]).
        // V is index 0, U is index 1.
        // If we iterate i from 0 to scc_cnt-1:
        // Process 0: Done.
        // Process 1: It has edge to 0. We need value of 0. 0 is ready.
        // So yes, iterating 0 to scc_cnt-1 allows pulling data from neighbors IF neighbors have lower indices.
        // But what if we have a back edge or forward edge relative to IDs?
        // In Tarjan, if U -> V and V is not visited, V gets lower ID.
        // If U -> V and V is visited (cross edge), V has lower ID.
        // So generally edges go High ID -> Low ID in Tarjan's DAG.
        // So we process 0, then 1... ensures children are ready.
        
        for (int next : dag_adj[i]) {
            // Note: In Tarjan's SCC DAG, edges usually go from Higher ID to Lower ID
            // BUT this is not guaranteed for all cases (e.g. disconnected components).
            // However, we need to update 'i' using 'next'. 
            // The logic above says 'next' should be processed before 'i'?
            // No, the loop direction is tricky.
            // Let's use the explicit topological order logic:
            // Tarjan's output order (0, 1, ... cnt-1) is Reverse Topological Order (Sink -> Source).
            // So if we iterate 0 to cnt-1, we are visiting Sinks first, then Sources.
            // When we visit Source U (index 10), it points to Sink V (index 0).
            // We want to update U using V.
            // V (index 0) was processed long ago.
            // So iterating 0 to cnt-1 works.
            
            // Wait, we iterate `i` (current SCC). We look at `next` (children).
            // Update `i` using `next`.
            scc_min[i] = min(scc_min[i], scc_min[next]);
            scc_max[i] = max(scc_max[i], scc_max[next]);
        }
    }

    // 6. Output for West nodes sorted by Y descending
    sort(west_points.begin(), west_points.end(), compareYDesc);

    for (const auto& wp : west_points) {
        int sid = scc_id[wp.id];
        if (scc_min[sid] > scc_max[sid]) {
            cout << "0\n";
        } else {
            cout << (scc_max[sid] - scc_min[sid] + 1) << "\n";
        }
    }

    return 0;
}