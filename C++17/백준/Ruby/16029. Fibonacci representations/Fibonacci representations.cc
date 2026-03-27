#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

const int MOD = 1000000007;
const int INF = 2000000007;

struct Matrix {
    int mat[2][2];
    Matrix() {
        mat[0][0] = 1; mat[0][1] = 0;
        mat[1][0] = 0; mat[1][1] = 1;
    }
    Matrix(int a, int b, int c, int d) {
        mat[0][0] = a % MOD; mat[0][1] = b % MOD;
        mat[1][0] = c % MOD; mat[1][1] = d % MOD;
    }
    Matrix operator*(const Matrix& o) const {
        Matrix res(0, 0, 0, 0);
        for(int i = 0; i < 2; ++i) {
            for(int j = 0; j < 2; ++j) {
                long long sum = 0;
                for(int k = 0; k < 2; ++k) {
                    sum += 1LL * mat[i][k] * o.mat[k][j];
                }
                res.mat[i][j] = sum % MOD;
            }
        }
        return res;
    }
};

struct Node {
    int min_val, max_val;
    Matrix M;
};

set<int> S;
struct Op {
    int type;
    int v;
};

vector<Op> ops;
vector<int> all_vals;
vector<Node> tree;

void add(int v) {
    if (S.count(v)) {
        S.erase(v);
        ops.push_back({0, v});
        if (v == 1) {
            add(2);
        } else if (v == 2) {
            add(3);
            add(1);
        } else {
            add(v + 1);
            add(v - 2);
        }
        return;
    }
    if (S.count(v - 1)) {
        S.erase(v - 1);
        ops.push_back({0, v - 1});
        add(v + 1);
        return;
    }
    if (S.count(v + 1)) {
        S.erase(v + 1);
        ops.push_back({0, v + 1});
        add(v + 2);
        return;
    }
    S.insert(v);
    ops.push_back({1, v});
    all_vals.push_back(v);
}

void build(int node, int l, int r) {
    tree[node].min_val = INF;
    tree[node].max_val = -INF;
    tree[node].M = Matrix();
    if (l == r) return;
    int mid = l + (r - l) / 2;
    build(node * 2, l, mid);
    build(node * 2 + 1, mid + 1, r);
}

void update(int node, int l, int r, int idx, int type, int val) {
    if (l == r) {
        if (type == 1) {
            tree[node].min_val = val;
            tree[node].max_val = val;
            tree[node].M = Matrix();
        } else {
            tree[node].min_val = INF;
            tree[node].max_val = -INF;
            tree[node].M = Matrix();
        }
        return;
    }
    int mid = l + (r - l) / 2;
    if (idx <= mid) update(node * 2, l, mid, idx, type, val);
    else update(node * 2 + 1, mid + 1, r, idx, type, val);

    int ls = node * 2, rs = node * 2 + 1;
    if (tree[ls].min_val == INF && tree[rs].min_val == INF) {
        tree[node].min_val = INF;
        tree[node].max_val = -INF;
        tree[node].M = Matrix();
    } else if (tree[ls].min_val == INF) {
        tree[node] = tree[rs];
    } else if (tree[rs].min_val == INF) {
        tree[node] = tree[ls];
    } else {
        tree[node].min_val = tree[ls].min_val;
        tree[node].max_val = tree[rs].max_val;
        int d = tree[rs].min_val - tree[ls].max_val;
        Matrix T(1, 1, (d - 1) / 2, d / 2);
        tree[node].M = tree[rs].M * T * tree[ls].M;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n;
    if (!(cin >> n)) return 0;
    
    ops.reserve(1000000);
    all_vals.reserve(1000000);
    
    for (int i = 0; i < n; i++) {
        int a;
        cin >> a;
        add(a);
        ops.push_back({2, 0});
    }

    sort(all_vals.begin(), all_vals.end());
    all_vals.erase(unique(all_vals.begin(), all_vals.end()), all_vals.end());

    int U = all_vals.size();
    if (U == 0) {
        for (int i = 0; i < n; i++) cout << 0 << "\n";
        return 0;
    }
    
    tree.resize(4 * U + 1);
    build(1, 1, U);

    for (auto& op : ops) {
        if (op.type == 2) {
            if (tree[1].min_val == INF) {
                cout << 0 << "\n";
            } else {
                int c1 = tree[1].min_val;
                Matrix T(1, 1, (c1 - 1) / 2, c1 / 2);
                Matrix M_final = tree[1].M * T;
                int ans = (M_final.mat[0][0] + M_final.mat[1][0]) % MOD;
                cout << ans << "\n";
            }
        } else {
            int idx = lower_bound(all_vals.begin(), all_vals.end(), op.v) - all_vals.begin() + 1;
            update(1, 1, U, idx, op.type, op.v);
        }
    }
    return 0;
}