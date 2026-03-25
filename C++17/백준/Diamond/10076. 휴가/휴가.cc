#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

const int MAXN = 100005;
long long count_seg[MAXN * 4];
long long sum_seg[MAXN * 4];

vector<int> compressed;
int M; 

void update_seg(int node, int l, int r, int pos, int cnt_diff, long long val_diff) {
    count_seg[node] += cnt_diff;
    sum_seg[node] += val_diff;
    if (l == r) return;
    int mid = (l + r) / 2;
    if (pos <= mid) update_seg(node * 2, l, mid, pos, cnt_diff, val_diff);
    else update_seg(node * 2 + 1, mid + 1, r, pos, cnt_diff, val_diff);
}

long long query(int node, int l, int r, int k) {
    if (k <= 0) return 0;
    if (count_seg[node] <= k) return sum_seg[node];
    if (l == r) return 1LL * k * compressed[l];
    int mid = (l + r) / 2;
    if (count_seg[node * 2 + 1] >= k) return query(node * 2 + 1, mid + 1, r, k);
    else return sum_seg[node * 2 + 1] + query(node * 2, l, mid, k - count_seg[node * 2 + 1]);
}

int curr_L, curr_R;

void add(int idx, const vector<int>& arr) {
    int v = arr[idx];
    int pos = lower_bound(compressed.begin(), compressed.end(), v) - compressed.begin();
    update_seg(1, 0, M - 1, pos, 1, v);
}

void remove(int idx, const vector<int>& arr) {
    int v = arr[idx];
    int pos = lower_bound(compressed.begin(), compressed.end(), v) - compressed.begin();
    update_seg(1, 0, M - 1, pos, -1, -v);
}

void move_pointers(int target_L, int target_R, const vector<int>& arr) {
    while (curr_L > target_L) { curr_L--; add(curr_L, arr); }
    while (curr_R < target_R) { curr_R++; add(curr_R, arr); }
    while (curr_L < target_L) { remove(curr_L, arr); curr_L++; }
    while (curr_R > target_R) { remove(curr_R, arr); curr_R--; }
}

long long ans_max = 0;

void solve(int y_min, int y_max, int x_min, int x_max, int start, int d, const vector<int>& arr) {
    if (y_min > y_max) return;
    int y_mid = (y_min + y_max) / 2;

    int best_x = x_min;
    long long best_val = -1;

    for (int x = x_min; x <= x_max; ++x) {
        move_pointers(start - x, start + y_mid, arr);
        int k = d - 2 * x - y_mid; // 남은 일수 (선택 가능한 관광지 수)
        
        long long current_val = 0;
        if (k > 0) current_val = query(1, 0, M - 1, k);
        
        if (current_val > best_val) {
            best_val = current_val;
            best_x = x;
        }
    }

    ans_max = max(ans_max, best_val);

    solve(y_min, y_mid - 1, best_x, x_max, start, d, arr);
    solve(y_mid + 1, y_max, x_min, best_x, start, d, arr);
}

void solve_half(int start, int d, const vector<int>& arr) {
    curr_L = start + 1;
    curr_R = start;
    memset(count_seg, 0, sizeof(count_seg));
    memset(sum_seg, 0, sizeof(sum_seg));
    
    solve(0, arr.size() - 1 - start, 0, start, start, d, arr);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, start, d;
    if (!(cin >> n >> start >> d)) return 0;

    vector<int> A(n);
    for (int i = 0; i < n; i++) {
        cin >> A[i];
        compressed.push_back(A[i]);
    }

    sort(compressed.begin(), compressed.end());
    compressed.erase(unique(compressed.begin(), compressed.end()), compressed.end());
    M = compressed.size();

    solve_half(start, d, A);

    vector<int> rev_A = A;
    reverse(rev_A.begin(), rev_A.end());
    int rev_start = n - 1 - start;

    solve_half(rev_start, d, rev_A);

    cout << ans_max << "\n";

    return 0;
}