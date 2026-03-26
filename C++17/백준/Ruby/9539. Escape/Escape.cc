#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const long long INF = 1e15;

// 우회 경로의 필요 체력과 순수 이득을 담는 구조체
struct Excursion {
    long long req, gain;
    bool operator<(const Excursion& other) const {
        if (req != other.req) return req > other.req; // 필요 체력이 작은 순서대로 우선
        return gain < other.gain; // 필요 체력이 같다면 이득이 큰 것을 우선
    }
};

void solve() {
    int n, t;
    cin >> n >> t;
    
    vector<long long> W(n + 1);
    for (int i = 1; i <= n; ++i) {
        cin >> W[i];
    }

    vector<vector<int>> adj(n + 1);
    for (int i = 0; i < n - 1; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // 트리 구조 파악을 위한 BFS 처리 (스택 오버플로우 방지)
    vector<int> order;
    vector<int> parent(n + 1, 0);
    queue<int> q;

    q.push(1);
    order.push_back(1);
    parent[1] = 0;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : adj[u]) {
            if (v != parent[u]) {
                parent[v] = u;
                q.push(v);
                order.push_back(v);
            }
        }
    }

    vector<priority_queue<Excursion>> pq(n + 1);

    // 단말 노드부터 역순으로 위로 올라가며 큐 병합
    for (int i = n - 1; i >= 0; --i) {
        int u = order[i];

        // 목표 지점에 도달하면 무한대의 보상 추가
        if (u == t) {
            pq[u].push({0, INF});
        }

        for (int v : adj[u]) {
            if (v != parent[u]) {
                // Small-to-Large Trick 병합
                if (pq[u].size() < pq[v].size()) {
                    swap(pq[u], pq[v]);
                }
                while (!pq[v].empty()) {
                    pq[u].push(pq[v].top());
                    pq[v].pop();
                }
            }
        }

        long long R = W[u] >= 0 ? 0 : -W[u];
        long long G = W[u];
        
        // 몬스터(음수 방)를 뚫기 위해 뒤의 보상들을 당겨오거나,
        // 위상적으로 몬스터 뒤에 종속되어야 하는 보상들을 강제로 묶음
        while (!pq[u].empty()) {
            Excursion top = pq[u].top();
            if (G <= 0 || top.req <= R) {
                pq[u].pop();
                R = max(R, top.req - G);
                G += top.gain;
            } else {
                break;
            }
        }
        
        // 최종적으로 이득이 있는 패키지만 부모에게 전달
        if (G > 0) {
            pq[u].push({R, G});
        }
    }

    // 1번 방(루트)에서 쌓여있는 선택지들을 그리디하게 챙겨봄
    long long current_HP = 0;
    while (!pq[1].empty()) {
        Excursion top = pq[1].top();
        pq[1].pop();
        
        if (current_HP < top.req) break;
        current_HP += top.gain;
    }

    // 최종적으로 무한대 보상 패키지를 무사히 열었는지 확인
    if (current_HP >= INF / 2) {
        cout << "escaped\n";
    } else {
        cout << "trapped\n";
    }
}

int main() {
    // 입출력 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        while (T--) {
            solve();
        }
    }
    return 0;
}