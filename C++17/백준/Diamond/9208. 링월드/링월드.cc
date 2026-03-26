#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Interval {
    long long x, y;
    bool operator<(const Interval& other) const {
        return y < other.y;
    }
};

int find_set(int v, vector<int>& parent_node) {
    if (v == parent_node[v])
        return v;
    return parent_node[v] = find_set(parent_node[v], parent_node);
}

void solve() {
    long long M;
    int N;
    if (!(cin >> M >> N)) return;
    
    vector<Interval> intervals;
    intervals.reserve(2 * N);
    
    bool impossible_early = false;
    if (N > M) {
        impossible_early = true;
    }
    
    for (int i = 0; i < N; ++i) {
        long long u, v;
        cin >> u >> v;
        if (impossible_early) continue;
        
        // 원형 구간을 두 배 길이의 직선으로 펼침
        long long len = (v - u + M) % M;
        long long y1 = u + len;
        intervals.push_back({u, y1});
        intervals.push_back({u + M, y1 + M});
    }
    
    if (impossible_early) {
        cout << "NO\n";
        return;
    }
    
    // 시작점 좌표들만 모아서 좌표 압축
    vector<long long> coords;
    coords.reserve(2 * N + 1);
    for (const auto& iv : intervals) {
        coords.push_back(iv.x);
    }
    coords.push_back(3000000005LL); // 최대 좌표 범위(3M)보다 큰 무한대 값 추가
    
    sort(coords.begin(), coords.end());
    coords.erase(unique(coords.begin(), coords.end()), coords.end());
    
    int K = coords.size() - 1;
    vector<int> parent_node(K + 1);
    vector<long long> capacity_val(K + 1);
    
    // 분리 집합 및 수용력 초기화
    for (int i = 0; i < K; ++i) {
        parent_node[i] = i;
        capacity_val[i] = coords[i + 1] - coords[i];
    }
    parent_node[K] = K; 
    capacity_val[K] = 0;
    
    // 구간을 끝점(y) 기준으로 오름차순 정렬
    sort(intervals.begin(), intervals.end());
    
    bool possible = true;
    for (const auto& iv : intervals) {
        int s = lower_bound(coords.begin(), coords.end(), iv.x) - coords.begin();
        int idx = find_set(s, parent_node);
        
        if (idx == K) {
            possible = false;
            break;
        }
        
        // 현재 블록에서 할당 가능한 가장 작은 값 계산
        long long P = coords[idx + 1] - capacity_val[idx];
        if (P > iv.y) {
            possible = false;
            break;
        }
        
        capacity_val[idx]--;
        // 수용력이 모두 소진되었다면 다음 블록으로 Union
        if (capacity_val[idx] == 0) {
            parent_node[idx] = find_set(idx + 1, parent_node);
        }
    }
    
    if (possible) cout << "YES\n";
    else cout << "NO\n";
}

int main() {
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