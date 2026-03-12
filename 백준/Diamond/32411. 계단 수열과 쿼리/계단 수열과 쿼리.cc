#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct Node {
    int pref[11];
    int suff[11];
    int max_len[11];
    int len;
    
    Node() {
        len = 0;
        for(int i = 1; i <= 10; ++i) pref[i] = suff[i] = max_len[i] = 0;
    }
};

Node tree[400005];
long long D[100005];

Node merge_node(const Node& L, const Node& R) {
    Node res;
    res.len = L.len + R.len;
    for(int k = 1; k <= 10; ++k) {
        res.pref[k] = L.pref[k] + (L.pref[k] == L.len ? R.pref[k] : 0);
        res.suff[k] = R.suff[k] + (R.suff[k] == R.len ? L.suff[k] : 0);
        res.max_len[k] = max({L.max_len[k], R.max_len[k], L.suff[k] + R.pref[k]});
    }
    return res;
}

void build(int node, int start, int end) {
    if (start == end) {
        tree[node].len = 1;
        for(int k = 1; k <= 10; ++k) {
            bool match = (abs(D[start]) == k);
            tree[node].pref[k] = tree[node].suff[k] = tree[node].max_len[k] = match ? 1 : 0;
        }
        return;
    }
    int mid = (start + end) / 2;
    build(2 * node, start, mid);
    build(2 * node + 1, mid + 1, end);
    tree[node] = merge_node(tree[2 * node], tree[2 * node + 1]);
}

void update(int node, int start, int end, int idx) {
    if (start == end) {
        for(int k = 1; k <= 10; ++k) {
            bool match = (abs(D[start]) == k);
            tree[node].pref[k] = tree[node].suff[k] = tree[node].max_len[k] = match ? 1 : 0;
        }
        return;
    }
    int mid = (start + end) / 2;
    if (idx <= mid) update(2 * node, start, mid, idx);
    else update(2 * node + 1, mid + 1, end, idx);
    tree[node] = merge_node(tree[2 * node], tree[2 * node + 1]);
}

Node query(int node, int start, int end, int l, int r) {
    if (l <= start && end <= r) return tree[node];
    int mid = (start + end) / 2;
    if (r <= mid) return query(2 * node, start, mid, l, r);
    if (l > mid) return query(2 * node + 1, mid + 1, end, l, r);
    return merge_node(query(2 * node, start, mid, l, r), query(2 * node + 1, mid + 1, end, l, r));
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n, q;
    if (!(cin >> n >> q)) return 0;
    
    vector<long long> a(n + 1);
    for(int i = 1; i <= n; ++i) {
        cin >> a[i];
    }
    
    if (n > 1) {
        for(int i = 1; i < n; ++i) {
            D[i] = a[i+1] - a[i];
        }
        build(1, 1, n - 1);
    }
    
    for(int i = 0; i < q; ++i) {
        int type;
        cin >> type;
        if (type == 1) {
            int l, r;
            long long x;
            cin >> l >> r >> x;
            if (l > 1) {
                D[l - 1] += x;
                update(1, 1, n - 1, l - 1);
            }
            if (r < n) {
                D[r] -= x;
                update(1, 1, n - 1, r);
            }
        } else if (type == 2) {
            int l, r, k;
            cin >> l >> r >> k;
            if (l == r) {
                cout << 1 << "\n";
            } else {
                Node res = query(1, 1, n - 1, l, r - 1);
                cout << res.max_len[k] + 1 << "\n";
            }
        }
    }
    
    return 0;
}