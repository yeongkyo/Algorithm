#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

// 최대 노드 수 (N + M <= 300,000)
const int MAX = 300005;

struct Edge {
    int to;
    long long weight;
};

vector<Edge> adj[MAX];
priority_queue<long long>* pq[MAX]; // 각 노드의 기울기 변화점들을 관리할 PQ
int leaf_cnt[MAX]; // 서브트리의 리프 노드 개수

void solve(int u) {
    // 1. 리프 노드인 경우 (자식이 없음)
    if (adj[u].empty()) {
        leaf_cnt[u] = 1;
        pq[u] = new priority_queue<long long>();
        pq[u]->push(0);
        pq[u]->push(0);
        return;
    }

    // 2. 내부 노드인 경우: 자식들의 PQ를 병합
    leaf_cnt[u] = 0;
    // 가장 큰 PQ를 찾아서 거기에 합치기 위한 포인터 초기화
    pq[u] = nullptr; 

    for (auto& edge : adj[u]) {
        int v = edge.to;
        long long w = edge.weight;
        
        solve(v); // 자식 노드 재귀 호출
        
        while (pq[v]->size() > leaf_cnt[v] + 1) {
            pq[v]->pop();
        }
        
        // 2) 간선 길이 w만큼 이동 (Shift logic)
        long long r = pq[v]->top(); pq[v]->pop();
        long long l = pq[v]->top(); pq[v]->pop();
        
        pq[v]->push(r + w);
        pq[v]->push(l + w);
        
        // 3) 병합 (Swap to larger)
        leaf_cnt[u] += leaf_cnt[v];
        if (pq[u] == nullptr) {
            pq[u] = pq[v];
        } else {
            if (pq[u]->size() < pq[v]->size()) {
                swap(pq[u], pq[v]);
            }
            // 작은 쪽(v)을 큰 쪽(u)에 다 털어넣음
            while (!pq[v]->empty()) {
                pq[u]->push(pq[v]->top());
                pq[v]->pop();
            }
            delete pq[v]; // 메모리 해제
        }
    }
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    long long total_weight_sum = 0;

    // 입력 처리
    // i는 2번 노드부터 N+M번 노드까지
    for (int i = 2; i <= n + m; i++) {
        int p;
        long long c;
        cin >> p >> c;
        adj[p].push_back({i, c});
        total_weight_sum += c;
    }

    // 루트(1번)에서 시작
    solve(1);

    int needed = leaf_cnt[1];
    while (pq[1]->size() > needed) {
        pq[1]->pop();
    }

    long long reduction = 0;
    while (!pq[1]->empty()) {
        reduction += pq[1]->top();
        pq[1]->pop();
    }

    // 최종 답: 전체 가중치 합 - (기울기 0 이하 구간의 변곡점 합)
    cout << total_weight_sum - reduction << "\n";

    return 0;
}