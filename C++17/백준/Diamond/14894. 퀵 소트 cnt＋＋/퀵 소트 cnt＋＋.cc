#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef long long ll;

const int MAX = 500005;
int n;
int a[MAX], pos[MAX];
int root[MAX];

struct Node {
    int count;
    int l, r;
} tree[MAX * 25]; 

int node_cnt = 0;

int update(int prev, int start, int end, int val) {
    int curr = ++node_cnt;
    tree[curr] = tree[prev];
    tree[curr].count++;
    if (start == end) return curr;

    int mid = (start + end) / 2;
    if (val <= mid) tree[curr].l = update(tree[prev].l, start, mid, val);
    else tree[curr].r = update(tree[prev].r, mid + 1, end, val);
    return curr;
}

int query(int nodeL, int nodeR, int start, int end, int k) {
    if (start == end) return start;
    int mid = (start + end) / 2;
    int diff = tree[tree[nodeR].l].count - tree[tree[nodeL].l].count;
    if (k <= diff) return query(tree[nodeL].l, tree[nodeR].l, start, mid, k);
    else return query(tree[nodeL].r, tree[nodeR].r, mid + 1, end, k - diff);
}

ll total_cnt = 0;

void solve(int valL, int valR) {
    if (valL >= valR) return;

    int size = valR - valL + 1;
    total_cnt += size;
    int k = (size - 1) / 2 + 1;
    int origin_idx = query(root[valL - 1], root[valR], 0, n - 1, k);
    int pivot_val = a[origin_idx];

    solve(valL, pivot_val - 1);
    solve(pivot_val + 1, valR);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    for (int i = 0; i < n; i++) {
        cin >> a[i];
        pos[a[i]] = i;
    }

    for (int i = 1; i <= n; i++) {
        root[i] = update(root[i - 1], 0, n - 1, pos[i]);
    }

    solve(1, n);

    cout << total_cnt << "\n";

    return 0;
}