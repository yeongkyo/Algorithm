#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef long long ll;

struct Edge {
    int to;
    int w_out; // u -> v 비용
    int w_in;  // v -> u 비용
};

int N;
vector<vector<Edge>> adj;
vector<ll> sz; // 서브트리 크기
vector<ll> dp; // 1번 노드 기준 서브트리 거리 합 (1차 DFS용)
ll min_cost = -1;

// 1차 DFS: 임의의 루트(1번)를 기준으로 서브트리 크기(sz)와 서브트리 내부 거리 합(dp) 계산
void dfs_bottom_up(int u, int p) {
    sz[u] = 1;
    dp[u] = 0;
    
    for (const auto& e : adj[u]) {
        if (e.to != p) {
            dfs_bottom_up(e.to, u);
            sz[u] += sz[e.to];
            // u에서 e.to 서브트리 내의 노드들로 가려면
            // u->e.to 간선을 'e.to 서브트리 노드 개수'만큼 지나가야 함
            dp[u] += dp[e.to] + (sz[e.to] * e.w_out);
        }
    }
}

// 2차 DFS: Rerooting을 통해 모든 노드를 루트로 했을 때의 총 비용 계산
void dfs_top_down(int u, int p, ll current_total_cost) {
    // 현재 노드 u를 우체국으로 했을 때의 비용으로 정답 갱신
    if (min_cost == -1 || current_total_cost < min_cost) {
        min_cost = current_total_cost;
    }

    for (const auto& e : adj[u]) {
        if (e.to != p) {
            // 루트를 u에서 e.to(자식)로 옮길 때 비용 변화 계산
            
            // 1. e.to의 서브트리에 있는 노드들은 u->e.to 거리만큼 가까워짐
            //    감소 비용 = (e.to 서브트리 크기) * (u->e.to 가중치)
            ll decrease = sz[e.to] * e.w_out;

            // 2. 그 외 나머지 노드들(전체 N - e.to 서브트리 크기)은 e.to->u 거리만큼 멀어짐
            //    증가 비용 = (N - e.to 서브트리 크기) * (e.to->u 가중치)
            //    주의: 도로는 방향성이 있으므로 돌아오는 비용(e.w_in)을 사용해야 함
            ll increase = (N - sz[e.to]) * e.w_in;

            ll next_total_cost = current_total_cost - decrease + increase;

            dfs_top_down(e.to, u, next_total_cost);
        }
    }
}

int main() {
    // 입출력 속도 향상
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N)) return 0;

    adj.resize(N + 1);
    sz.resize(N + 1);
    dp.resize(N + 1);

    for (int i = 0; i < N - 1; ++i) {
        int u, v, w_u, w_v;
        cin >> u >> v >> w_u >> w_v;
        // u -> v 비용 w_u, 돌아오는 비용 w_v
        adj[u].push_back({v, w_u, w_v});
        // v -> u 비용 w_v, 돌아오는 비용 w_u
        adj[v].push_back({u, w_v, w_u});
    }

    // 1. Bottom-Up으로 1번 노드 기준 정보 계산
    dfs_bottom_up(1, 0);

    // 2. Top-Down으로 루트를 옮겨가며 최솟값 찾기
    // 초기값은 1번 노드가 루트일 때의 dp값
    dfs_top_down(1, 0, dp[1]);

    cout << min_cost << "\n";

    return 0;
}