#include <iostream>

using namespace std;

const int MOD = 1000000007;

struct Matrix {
    int mat[2][2];
    Matrix() {
        mat[0][0] = 1; mat[0][1] = 0;
        mat[1][0] = 0; mat[1][1] = 1;
    }
    Matrix(int a, int b, int c, int d) {
        mat[0][0] = a; mat[0][1] = b;
        mat[1][0] = c; mat[1][1] = d;
    }
    Matrix operator*(const Matrix& o) const {
        return Matrix(
            (1LL * mat[0][0] * o.mat[0][0] + 1LL * mat[0][1] * o.mat[1][0]) % MOD,
            (1LL * mat[0][0] * o.mat[0][1] + 1LL * mat[0][1] * o.mat[1][1]) % MOD,
            (1LL * mat[1][0] * o.mat[0][0] + 1LL * mat[1][1] * o.mat[1][0]) % MOD,
            (1LL * mat[1][0] * o.mat[0][1] + 1LL * mat[1][1] * o.mat[1][1]) % MOD
        );
    }
};

unsigned int xorshift() {
    static unsigned int x = 123456789;
    static unsigned int y = 362436069;
    static unsigned int z = 521288629;
    static unsigned int w = 88675123;
    unsigned int t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

struct Node {
    int key;
    unsigned int priority;
    int min_val, max_val;
    Matrix M;
    int l, r;
};

const int MAX_NODES = 5000005;
Node pool[MAX_NODES];
int node_cnt = 0;
int trash[MAX_NODES];
int trash_cnt = 0;

int allocNode(int key) {
    int id = trash_cnt > 0 ? trash[--trash_cnt] : ++node_cnt;
    pool[id].key = key;
    pool[id].priority = xorshift();
    pool[id].min_val = pool[id].max_val = key;
    pool[id].M = Matrix();
    pool[id].l = pool[id].r = 0;
    return id;
}

void freeNode(int id) {
    if (id) trash[trash_cnt++] = id;
}

void update(int t) {
    if (!t) return;
    int l = pool[t].l;
    int r = pool[t].r;

    pool[t].min_val = l ? pool[l].min_val : pool[t].key;
    pool[t].max_val = r ? pool[r].max_val : pool[t].key;

    Matrix M(1, 0, 0, 1);
    if (l) {
        int d = pool[t].key - pool[l].max_val;
        Matrix T(1, 1, (d - 1) / 2, d / 2);
        M = T * pool[l].M;
    }
    if (r) {
        int d = pool[r].min_val - pool[t].key;
        Matrix T(1, 1, (d - 1) / 2, d / 2);
        M = pool[r].M * T * M;
    }
    pool[t].M = M;
}

void split(int t, int key, int &l, int &r) {
    if (!t) { l = r = 0; return; }
    if (pool[t].key < key) {
        split(pool[t].r, key, pool[t].r, r);
        l = t;
    } else {
        split(pool[t].l, key, l, pool[t].l);
        r = t;
    }
    update(t);
}

void merge(int &t, int l, int r) {
    if (!l || !r) { t = l ? l : r; return; }
    if (pool[l].priority > pool[r].priority) {
        merge(pool[l].r, pool[l].r, r);
        t = l;
    } else {
        merge(pool[r].l, l, pool[r].l);
        t = r;
    }
    update(t);
}

bool find(int t, int key) {
    while (t) {
        if (pool[t].key == key) return true;
        if (pool[t].key < key) t = pool[t].r;
        else t = pool[t].l;
    }
    return false;
}

void insert(int &t, int key) {
    int l, r;
    split(t, key, l, r);
    int mid = allocNode(key);
    merge(t, l, mid);
    merge(t, t, r);
}

void erase(int &t, int key) {
    int l, mid, r;
    split(t, key, l, r);
    split(r, key + 1, mid, r);
    freeNode(mid);
    merge(t, l, r);
}

int root = 0;

void add(int v) {
    if (find(root, v)) {
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
    if (find(root, v - 1)) {
        erase(root, v - 1);
        add(v + 1);
        return;
    }
    if (find(root, v + 1)) {
        erase(root, v + 1);
        add(v + 2);
        return;
    }
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
            int c1 = pool[root].min_val;
            Matrix T(1, 1, (c1 - 1) / 2, c1 / 2);
            Matrix M_final = pool[root].M * T;
            int ans = (M_final.mat[0][0] + M_final.mat[1][0]) % MOD;
            cout << ans << "\n";
        }
    }
    return 0;
}