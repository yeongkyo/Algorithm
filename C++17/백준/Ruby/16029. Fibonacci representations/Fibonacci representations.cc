#include <iostream>
#include <set>
#include <random>

using namespace std;

const int MOD = 1000000007;

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

mt19937 rng(1337);

struct Node {
    int key;
    int priority;
    int min_val, max_val;
    Matrix M;
    Node *l, *r;
    Node(int k) : key(k), priority(rng()), min_val(k), max_val(k), M(), l(NULL), r(NULL) {}
};

typedef Node* pNode;

void update(pNode t) {
    if (!t) return;
    t->min_val = t->key;
    t->max_val = t->key;
    t->M = Matrix();

    if (t->l) {
        int d = t->key - t->l->max_val;
        Matrix T(1, 1, (d - 1) / 2, d / 2);
        t->M = t->M * T * t->l->M;
        t->min_val = t->l->min_val;
    }
    if (t->r) {
        int d = t->r->min_val - t->max_val;
        Matrix T(1, 1, (d - 1) / 2, d / 2);
        t->M = t->r->M * T * t->M;
        t->max_val = t->r->max_val;
    }
}

void split(pNode t, int key, pNode &l, pNode &r) {
    if (!t) { l = r = NULL; return; }
    if (t->key < key) {
        split(t->r, key, t->r, r);
        l = t;
    } else {
        split(t->l, key, l, t->l);
        r = t;
    }
    update(t);
}

void merge(pNode &t, pNode l, pNode r) {
    if (!l || !r) { t = l ? l : r; return; }
    if (l->priority > r->priority) {
        merge(l->r, l->r, r);
        t = l;
    } else {
        merge(r->l, l, r->l);
        t = r;
    }
    update(t);
}

void insert(pNode &t, int key) {
    pNode l, r;
    split(t, key, l, r);
    pNode mid = new Node(key);
    merge(t, l, mid);
    merge(t, t, r);
}

void erase(pNode &t, int key) {
    pNode l, mid, r;
    split(t, key, l, r);
    split(r, key + 1, mid, r);
    if (mid) delete mid;
    merge(t, l, r);
}

pNode root = NULL;
set<int> S;

void add(int v) {
    if (S.count(v)) {
        S.erase(v);
        erase(root, v);
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
        erase(root, v - 1);
        add(v + 1);
        return;
    }
    if (S.count(v + 1)) {
        S.erase(v + 1);
        erase(root, v + 1);
        add(v + 2);
        return;
    }
    S.insert(v);
    insert(root, v);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    for (int i = 0; i < n; i++) {
        int a;
        cin >> a;
        add(a);

        if (!root) {
            cout << 0 << "\n";
        } else {
            int c1 = root->min_val;
            Matrix T(1, 1, (c1 - 1) / 2, c1 / 2);
            Matrix M_final = root->M * T;
            int ans = (M_final.mat[0][0] + M_final.mat[1][0]) % MOD;
            cout << ans << "\n";
        }
    }
    return 0;
}