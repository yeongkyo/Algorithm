#include <iostream>
#include <vector>
#include <random>

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
        res.mat[0][0] = (mat[0][0] * other.mat[0][0] + mat[0][1] * other.mat[1][0]) % MOD;
        res.mat[0][1] = (mat[0][0] * other.mat[0][1] + mat[0][1] * other.mat[1][1]) % MOD;
        res.mat[1][0] = (mat[1][0] * other.mat[0][0] + mat[1][1] * other.mat[1][0]) % MOD;
        res.mat[1][1] = (mat[1][0] * other.mat[0][1] + mat[1][1] * other.mat[1][1]) % MOD;
        return res;
    }
};

struct Node {
    int val;
    unsigned int pri;
    int ls, rs;
    int min_val, max_val;
    Matrix mat;
} tree[400005];

int root = 0, node_cnt = 0;
int free_list[400005], free_cnt = 0;
mt19937 rng(1337);

int new_node(int val) {
    int u = free_cnt > 0 ? free_list[--free_cnt] : ++node_cnt;
    tree[u].val = val;
    tree[u].pri = rng();
    tree[u].ls = tree[u].rs = 0;
    tree[u].min_val = tree[u].max_val = val;
    tree[u].mat = Matrix();
    return u;
}

void push_up(int u) {
    int ls = tree[u].ls;
    int rs = tree[u].rs;
    
    tree[u].min_val = ls ? tree[ls].min_val : tree[u].val;
    tree[u].max_val = rs ? tree[rs].max_val : tree[u].val;
    
    Matrix res; 
    if (ls) {
        res = tree[ls].mat;
        res = Matrix(tree[u].val - tree[ls].max_val) * res;
    }
    if (rs) {
        Matrix right_part = tree[rs].mat;
        right_part = right_part * Matrix(tree[rs].min_val - tree[u].val);
        res = right_part * res;
    }
    tree[u].mat = res;
}

void split(int u, int val, int& x, int& y) {
    if (!u) {
        x = y = 0;
        return;
    }
    if (tree[u].val <= val) {
        x = u;
        split(tree[u].rs, val, tree[u].rs, y);
    } else {
        y = u;
        split(tree[u].ls, val, x, tree[u].ls);
    }
    push_up(u);
}

int merge(int x, int y) {
    if (!x || !y) return x ? x : y;
    if (tree[x].pri > tree[y].pri) {
        tree[x].rs = merge(tree[x].rs, y);
        push_up(x);
        return x;
    } else {
        tree[y].ls = merge(x, tree[y].ls);
        push_up(y);
        return y;
    }
}

void insert_treap(int val) {
    int x, y;
    split(root, val, x, y);
    root = merge(merge(x, new_node(val)), y);
}

void collect_free(int u) {
    if (!u) return;
    free_list[free_cnt++] = u;
    collect_free(tree[u].ls);
    collect_free(tree[u].rs);
}

void erase_treap(int val) {
    int x, y, z;
    split(root, val - 1, x, y);
    split(y, val, y, z);
    collect_free(y);
    root = merge(x, z);
}

bool contains(int val) {
    int u = root;
    while (u) {
        if (tree[u].val == val) return true;
        if (val < tree[u].val) u = tree[u].ls;
        else u = tree[u].rs;
    }
    return false;
}

void add(int x) {
    if (x <= 0) return;
    if (contains(x)) {
        erase_treap(x);
        add(x + 1);
        if (x == 1) ;
        else if (x == 2) add(1);
        else add(x - 2);
        return;
    }
    if (contains(x + 1)) {
        erase_treap(x + 1);
        add(x + 2);
        return;
    }
    if (contains(x - 1)) {
        erase_treap(x - 1);
        add(x + 1);
        return;
    }
    insert_treap(x);
}

long long get_ans() {
    if (!root) return 0;
    Matrix total = tree[root].mat;
    Matrix first_gap(tree[root].min_val);
    total = total * first_gap;
    long long ans = (total.mat[0][0] + total.mat[1][0]) % MOD;
    return ans;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n;
    if (!(cin >> n)) return 0;
    for (int i = 0; i < n; ++i) {
        int a;
        cin >> a;
        add(a);
        cout << get_ans() << "\n";
    }
    return 0;
}