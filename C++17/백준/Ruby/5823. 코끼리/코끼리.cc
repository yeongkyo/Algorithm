#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int B = 400; // 블록(버킷)의 크기 (약 sqrt(150000))
int N, L, M;
vector<int> P; // 각 코끼리의 현재 위치

// 각 블록을 관리하는 구조체
struct Bucket {
    vector<int> p;     // 블록 내 코끼리들의 위치 (항상 정렬된 상태 유지)
    vector<int> id;    // 코끼리의 고유 번호
    vector<int> cams;  // i번째부터 블록 끝까지 덮는 데 필요한 카메라 수
    vector<int> reach; // i번째부터 덮었을 때, 마지막 카메라가 도달하는 최대 좌표

    // 블록 내부의 cams와 reach 배열을 O(K)에 재계산하는 함수
    void calc() {
        int sz = p.size();
        cams.resize(sz);
        reach.resize(sz);
        int j = sz - 1;
        
        // 뒤에서부터 역순으로 DP처럼 계산
        for (int i = sz - 1; i >= 0; --i) {
            int target = p[i] + L;
            // 현재 카메라(p[i] + L) 범위를 벗어나는 첫 번째 코끼리 j를 찾음
            while (j >= 0 && p[j] > target) j--;
            int next_idx = j + 1;
            
            if (next_idx < sz) {
                cams[i] = 1 + cams[next_idx];
                reach[i] = reach[next_idx];
            } else {
                cams[i] = 1;
                reach[i] = target;
            }
        }
    }
};

vector<Bucket> buckets;

// 전체 블록의 균형이 무너지는 것을 막기 위해 주기적으로 재배열하는 함수
void rebuild() {
    vector<pair<int, int>> all_e;
    all_e.reserve(N);
    for (const auto& b : buckets) {
        for (size_t i = 0; i < b.p.size(); ++i) {
            all_e.push_back({b.p[i], b.id[i]});
        }
    }
    
    buckets.clear();
    int num_buckets = (N + B - 1) / B;
    if (num_buckets == 0) num_buckets = 1; 
    buckets.resize(num_buckets);
    
    for (size_t i = 0; i < all_e.size(); ++i) {
        int b_idx = i / B;
        buckets[b_idx].p.push_back(all_e[i].first);
        buckets[b_idx].id.push_back(all_e[i].second);
    }
    
    for (auto& b : buckets) {
        if (!b.p.empty()) b.calc();
    }
}

// 전체 코끼리를 덮는 최소 카메라 수를 구하는 함수
int query() {
    int ans = 0;
    int last_reach = -1;
    
    for (const auto& b : buckets) {
        // 블록이 비어있거나, 이전 카메라의 도달 거리 안에 블록 전체가 포함되면 패스
        if (b.p.empty() || b.p.back() <= last_reach) continue;
        
        // 이전 카메라 범위를 벗어나는 첫 번째 코끼리를 이분 탐색으로 찾음
        auto it = upper_bound(b.p.begin(), b.p.end(), last_reach);
        int idx = distance(b.p.begin(), it);
        
        ans += b.cams[idx];
        last_reach = b.reach[idx];
    }
    return ans;
}

// 코끼리 이동 업데이트
void update(int id, int y) {
    int old_p = P[id];
    P[id] = y;

    // 1. 기존 위치에서 코끼리 삭제
    for (auto& b : buckets) {
        if (b.p.empty()) continue;
        if (old_p >= b.p.front() && old_p <= b.p.back()) {
            bool found = false;
            for (size_t i = 0; i < b.p.size(); ++i) {
                if (b.id[i] == id) {
                    b.p.erase(b.p.begin() + i);
                    b.id.erase(b.id.begin() + i);
                    found = true;
                    break;
                }
            }
            if (found) {
                b.calc(); // 삭제가 일어난 블록 재계산
                break;
            }
        }
    }

    // 2. 새로운 위치에 코끼리 삽입
    int target_b = buckets.size() - 1;
    for (int i = 0; i < (int)buckets.size(); ++i) {
        if (!buckets[i].p.empty() && y <= buckets[i].p.back()) {
            target_b = i;
            break;
        }
    }
    
    // 이분 탐색으로 적절한 위치를 찾아 삽입하여 정렬 상태 유지
    auto it = lower_bound(buckets[target_b].p.begin(), buckets[target_b].p.end(), y);
    int idx = distance(buckets[target_b].p.begin(), it);
    buckets[target_b].p.insert(buckets[target_b].p.begin() + idx, y);
    buckets[target_b].id.insert(buckets[target_b].id.begin() + idx, id);
    buckets[target_b].calc(); // 삽입이 일어난 블록 재계산
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> L >> M)) return 0;

    P.resize(N);
    for (int i = 0; i < N; ++i) {
        cin >> P[i];
    }

    // 초기 버킷 세팅
    int num_buckets = (N + B - 1) / B;
    if (num_buckets == 0) num_buckets = 1;
    buckets.resize(num_buckets);
    
    for (int i = 0; i < N; ++i) {
        int b_idx = i / B;
        buckets[b_idx].p.push_back(P[i]);
        buckets[b_idx].id.push_back(i);
    }
    for (auto& b : buckets) {
        if (!b.p.empty()) b.calc();
    }

    int update_cnt = 0;
    for (int i = 0; i < M; ++i) {
        int id, y;
        cin >> id >> y;
        
        update(id, y);
        update_cnt++;
        
        // 블록의 불균형을 막기 위해 주기적으로 전체 리빌딩
        if (update_cnt >= B) {
            rebuild();
            update_cnt = 0;
        }
        
        // 정답 출력
        cout << query() << "\n";
    }

    return 0;
}