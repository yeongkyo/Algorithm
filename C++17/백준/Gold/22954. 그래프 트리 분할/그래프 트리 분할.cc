#include <iostream>
#include <vector>

using namespace std;

struct Edge {
    int u, v, idx;
};

int N, M;
vector<pair<int, int>> adj[100005];
bool visited[100005];

vector<int> comp_nodes[3];
vector<Edge> comp_edges[3];

// 연결 요소를 찾으면서 스패닝 트리를 구성하는 DFS
void dfs_comp(int u, int c) {
    comp_nodes[c].push_back(u);
    visited[u] = true;
    for (auto& edge : adj[u]) {
        int v = edge.first;
        int idx = edge.second;
        if (!visited[v]) {
            comp_edges[c].push_back({u, v, idx});
            dfs_comp(v, c);
        }
    }
}

vector<pair<int, int>> tree_adj[100005];
int sub_size[100005];
int parent_edge[100005];

// 트리의 서브트리 크기를 구하는 DFS
void dfs_size(int u, int p) {
    sub_size[u] = 1;
    for (auto& edge : tree_adj[u]) {
        int v = edge.first;
        int idx = edge.second;
        if (v != p) {
            parent_edge[v] = idx;
            dfs_size(v, u);
            sub_size[u] += sub_size[v];
        }
    }
}

vector<int> t_nodes, t_edges;

// 잘라낸 간선을 제외하고 트리를 탐색해 정점과 간선을 수집
void get_tree(int u, int p, int cut_idx) {
    t_nodes.push_back(u);
    for (auto& edge : tree_adj[u]) {
        int v = edge.first;
        int idx = edge.second;
        if (v != p && idx != cut_idx) {
            t_edges.push_back(idx);
            get_tree(v, u, cut_idx);
        }
    }
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M)) return 0;

    for (int i = 1; i <= M; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back({v, i});
        adj[v].push_back({u, i});
    }

    // 정점이 2개 이하면 크기가 다른 두 트리로 나눌 수 없음 (1+1=2 이므로 크기 같음)
    if (N <= 2) {
        cout << -1 << "\n";
        return 0;
    }

    int C = 0;
    for (int i = 1; i <= N; ++i) {
        if (!visited[i]) {
            if (C == 2) {
                // 연결 요소가 3개 이상이면 분할 불가
                cout << -1 << "\n";
                return 0;
            }
            dfs_comp(i, C);
            C++;
        }
    }

    // 연결 요소가 2개인 경우
    if (C == 2) {
        if (comp_nodes[0].size() == comp_nodes[1].size()) {
            cout << -1 << "\n";
        } else {
            cout << comp_nodes[0].size() << " " << comp_nodes[1].size() << "\n";
            for (int u : comp_nodes[0]) cout << u << " "; cout << "\n";
            if (comp_edges[0].empty()) cout << "\n";
            else { for (auto& e : comp_edges[0]) cout << e.idx << " "; cout << "\n"; }

            for (int u : comp_nodes[1]) cout << u << " "; cout << "\n";
            if (comp_edges[1].empty()) cout << "\n";
            else { for (auto& e : comp_edges[1]) cout << e.idx << " "; cout << "\n"; }
        }
        return 0;
    }

    // 연결 요소가 1개인 경우 (C == 1)
    for (auto& e : comp_edges[0]) {
        tree_adj[e.u].push_back({e.v, e.idx});
        tree_adj[e.v].push_back({e.u, e.idx});
    }

    // 서브트리 크기 계산
    dfs_size(1, 0);

    int cut_v = -1;
    for (int i = 2; i <= N; ++i) {
        // 서브트리의 크기가 전체의 절반이 아닌 경우를 찾음
        if (sub_size[i] * 2 != N) {
            cut_v = i;
            break;
        }
    }

    // 조건을 만족하는 자를 간선이 없다면
    if (cut_v == -1) {
        cout << -1 << "\n";
        return 0;
    }

    int cut_idx = parent_edge[cut_v];

    // 첫 번째 트리 분리
    t_nodes.clear(); t_edges.clear();
    get_tree(cut_v, 0, cut_idx);
    vector<int> n1 = t_nodes;
    vector<int> e1 = t_edges;

    // 두 번째 트리 분리 (루트인 1번 정점부터 탐색)
    t_nodes.clear(); t_edges.clear();
    get_tree(1, 0, cut_idx);
    vector<int> n2 = t_nodes;
    vector<int> e2 = t_edges;

    // 결과 출력
    cout << n1.size() << " " << n2.size() << "\n";
    
    for (int u : n1) cout << u << " "; cout << "\n";
    if (e1.empty()) cout << "\n";
    else { for (int idx : e1) cout << idx << " "; cout << "\n"; }

    for (int u : n2) cout << u << " "; cout << "\n";
    if (e2.empty()) cout << "\n";
    else { for (int idx : e2) cout << idx << " "; cout << "\n"; }

    return 0;
}