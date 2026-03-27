#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

const int MOD = 1e9 + 7;

struct Matrix {
    long long mat[2][2];
    Matrix() {
        mat[0][0] = 1; mat[0][1] = 0;
        mat[1][0] = 0; mat[1][1] = 1;
    }
    Matrix(long long d) {
        mat[0][0] = 1;
        mat[0][1] = 1;
        mat[1][0] = ((d - 1) / 2) % MOD;
        mat[1][1] = (d / 2) % MOD;
    }
    Matrix operator*(const Matrix& other) const {
        Matrix res;
        long long a00 = mat[0][0], a01 = mat[0][1];
        long long a10 = mat[1][0], a11 = mat[1][1];
        long long b00 = other.mat[0][0], b01 = other.mat[0][1];
        long long b10 = other.mat[1][0], b11 = other.mat[1][1];
        
        res.mat[0][0] = (a00 * b00 + a01 * b10) % MOD;
        res.mat[0][1] = (a00 * b01 + a01 * b11) % MOD;
        res.mat[1][0] = (a10 * b00 + a11 * b10) % MOD;
        res.mat[1][1] = (a10 * b01 + a11 * b11) % MOD;
        return res;
    }
};

struct Node {
    bool active;
    int min_val, max_val;
    Matrix mat;
};

vector<Node> tree;
vector<int> vals;

void push_up(int u) {
    int lc = u << 1, rc = u << 1 | 1;
    if (!tree[lc].active && !tree[rc].active) {
        tree[u].active = false;
        return;
    }
    if (!tree[lc].active) {
        tree[u] = tree[rc];
        return;
    }
    if (!tree[rc].active) {
        tree[u] = tree[lc];
        return;
    }
    
    tree[u].active = true;
    tree[u].min_val = tree[lc].min_val;
    tree[u].max_val = tree[rc].max_val;
    
    int gap = tree[rc].min_val - tree[lc].max_val;
    tree[u].mat = tree[rc].mat * Matrix(gap) * tree[lc].mat;
}

void build(int u, int l, int r) {
    tree[u].active = false;
    if (l == r) {
        tree[u].min_val = tree[u].max_val = vals[l];
        tree[u].mat = Matrix();
        return;
    }
    int mid = l + (r - l) / 2;
    build(u << 1, l, mid);
    build(u << 1 | 1, mid + 1, r);
}

void update(int u, int l, int r, int idx, bool val) {
    if (l == r) {
        tree[u].active = val;
        return;
    }
    int mid = l + (r - l) / 2;
    if (idx <= mid) update(u << 1, l, mid, idx, val);
    else update(u << 1 | 1, mid + 1, r, idx, val);
    push_up(u);
}

long long get_ans() {
    if (!tree[1].active) return 0;
    Matrix total = tree[1].mat;
    Matrix first_gap(tree[1].min_val);
    total = total * first_gap;
    return (total.mat[0][0] + total.mat[1][0]) % MOD;
}

void sim_add(int x, set<int>& S, int current_query, vector<vector<pair<int, int>>>& ops_per_query) {
    if (x <= 0) return;
    if (S.count(x)) {
        S.erase(x);
        ops_per_query[current_query].push_back({x, -1});
        sim_add(x + 1, S, current_query, ops_per_query);
        if (x == 1) {}
        else if (x == 2) sim_add(1, S, current_query, ops_per_query);
        else sim_add(x - 2, S, current_query, ops_per_query);
        return;
    }
    if (S.count(x + 1)) {
        S.erase(x + 1);
        ops_per_query[current_query].push_back({x + 1, -1});
        sim_add(x + 2, S, current_query, ops_per_query);
        return;
    }
    if (S.count(x - 1)) {
        S.erase(x - 1);
        ops_per_query[current_query].push_back({x - 1, -1});
        sim_add(x + 1, S, current_query, ops_per_query);
        return;
    }
    S.insert(x);
    ops_per_query[current_query].push_back({x, 1});
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n;
    if (!(cin >> n)) return 0;
    
    vector<vector<pair<int, int>>> ops_per_query(n);
    set<int> S;
    
    for (int i = 0; i < n; ++i) {
        int a;
        cin >> a;
        sim_add(a, S, i, ops_per_query);
    }
    
    for (int i = 0; i < n; ++i) {
        for (auto& op : ops_per_query[i]) {
            vals.push_back(op.first);
        }
    }
    
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    
    int M = vals.size();
    if (M == 0) return 0;
    
    tree.resize(4 * M);
    build(1, 0, M - 1);
    
    for (int i = 0; i < n; ++i) {
        for (auto& op : ops_per_query[i]) {
            int idx = lower_bound(vals.begin(), vals.end(), op.first) - vals.begin();
            update(1, 0, M - 1, idx, op.second == 1);
        }
        cout << get_ans() << "\n";
    }
    
    return 0;
}