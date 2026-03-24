#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

const int MOD = 998244353;

struct Node {
    int max_val;
    int lazy;
};

int N;
vector<Node> tree;

void push_down(int node) {
    if (tree[node].lazy != -1) {
        tree[2 * node].max_val = max(tree[2 * node].max_val, tree[node].lazy);
        tree[2 * node].lazy = max(tree[2 * node].lazy, tree[node].lazy);
        tree[2 * node + 1].max_val = max(tree[2 * node + 1].max_val, tree[node].lazy);
        tree[2 * node + 1].lazy = max(tree[2 * node + 1].lazy, tree[node].lazy);
        tree[node].lazy = -1;
    }
}

void update_point(int node, int start, int end, int idx, int val) {
    if (start == end) {
        tree[node].max_val = val;
        tree[node].lazy = -1;
        return;
    }
    push_down(node);
    int mid = (start + end) / 2;
    if (idx <= mid) update_point(2 * node, start, mid, idx, val);
    else update_point(2 * node + 1, mid + 1, end, idx, val);
    tree[node].max_val = max(tree[2 * node].max_val, tree[2 * node + 1].max_val);
}

void update_chmax(int node, int start, int end, int l, int r, int val) {
    if (l > end || r < start) return;
    if (l <= start && end <= r) {
        tree[node].max_val = max(tree[node].max_val, val);
        tree[node].lazy = max(tree[node].lazy, val);
        return;
    }
    push_down(node);
    int mid = (start + end) / 2;
    update_chmax(2 * node, start, mid, l, r, val);
    update_chmax(2 * node + 1, mid + 1, end, l, r, val);
    tree[node].max_val = max(tree[2 * node].max_val, tree[2 * node + 1].max_val);
}

int query(int node, int start, int end, int idx) {
    if (start == end) return tree[node].max_val;
    push_down(node);
    int mid = (start + end) / 2;
    if (idx <= mid) return query(2 * node, start, mid, idx);
    else return query(2 * node + 1, mid + 1, end, idx);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N)) return 0;

    vector<int> A(N);
    for (int i = 0; i < N; i++) {
        cin >> A[i];
    }

    tree.assign(4 * N + 4, {-1, -1});
    vector<long long> V(N + 1, 0);
    vector<long long> W(N + 1, 0);
    
    set<int> valid_P;
    set<int> trapped_P;

    int x = A[0];
    V[x] = 1;
    long long SumV = 1;
    valid_P.insert(x);

    for (int i = 1; i < N; i++) {
        x = A[i];
        long long new_Vx = SumV;

        auto it = trapped_P.upper_bound(x);
        while (it != trapped_P.end()) {
            int P = *it;
            int M = query(1, 0, N - 1, P);
            V[M] = (V[M] + W[P]) % MOD;
            SumV = (SumV + W[P]) % MOD;
            valid_P.insert(M);
            W[P] = 0;
            it = trapped_P.erase(it);
        }

        auto it2 = valid_P.begin();
        while (it2 != valid_P.end() && *it2 < x) {
            int P = *it2;
            W[P] = V[P];
            V[P] = 0;
            SumV = (SumV - W[P] % MOD + MOD) % MOD;
            update_point(1, 0, N - 1, P, x);
            trapped_P.insert(P);
            it2 = valid_P.erase(it2);
        }

        if (!trapped_P.empty()) {
            update_chmax(1, 0, N - 1, 0, x, x);
        }

        V[x] = (V[x] + new_Vx) % MOD;
        SumV = (SumV + new_Vx) % MOD;
        valid_P.insert(x);
    }

    cout << SumV << "\n";

    return 0;
}