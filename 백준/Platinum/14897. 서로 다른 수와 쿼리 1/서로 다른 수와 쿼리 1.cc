#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MAX_N = 1000005;

int N, Q;
int A[MAX_N];
int last_pos[MAX_N]; // 압축된 좌표의 마지막 등장 위치
int bit[MAX_N];
int ans[MAX_N];
vector<pair<int, int>> queries[MAX_N]; // queries[r] = {l, query_index}
vector<int> compressed_val;

void update(int idx, int val) {
    for (; idx <= N; idx += idx & -idx)
        bit[idx] += val;
}

int query(int idx) {
    int sum = 0;
    for (; idx > 0; idx -= idx & -idx)
        sum += bit[idx];
    return sum;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> N;
    for (int i = 1; i <= N; i++) {
        cin >> A[i];
        compressed_val.push_back(A[i]);
    }

    // 좌표 압축
    sort(compressed_val.begin(), compressed_val.end());
    compressed_val.erase(unique(compressed_val.begin(), compressed_val.end()), compressed_val.end());

    for (int i = 1; i <= N; i++) {
        A[i] = lower_bound(compressed_val.begin(), compressed_val.end(), A[i]) - compressed_val.begin() + 1;
    }

    cin >> Q;
    for (int i = 0; i < Q; i++) {
        int l, r;
        cin >> l >> r;
        queries[r].push_back({l, i});
    }

    // 스위핑
    for (int i = 1; i <= N; i++) {
        int val = A[i];
        
        // 이전에 등장했던 수라면 이전 위치의 카운트를 제거
        if (last_pos[val] != 0) {
            update(last_pos[val], -1);
        }
        
        // 현재 위치에 카운트 추가
        update(i, 1);
        last_pos[val] = i;

        // 현재 위치(r)로 끝나는 쿼리 처리
        for (auto& q : queries[i]) {
            int l = q.first;
            int idx = q.second;
            ans[idx] = query(i) - query(l - 1);
        }
    }

    for (int i = 0; i < Q; i++) {
        cout << ans[i] << "\n";
    }

    return 0;
}