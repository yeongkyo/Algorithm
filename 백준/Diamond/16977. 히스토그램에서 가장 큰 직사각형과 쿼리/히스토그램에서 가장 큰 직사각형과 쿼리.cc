#include <bits/stdc++.h>
using namespace std;

struct FastScanner {
    static inline int gc() { return getchar_unlocked(); }

    template <typename T>
    bool readInt(T &out) {
        int c = gc();
        while (c != EOF && c <= ' ') c = gc();
        if (c == EOF) return false;

        T sign = 1;
        if (c == '-') { sign = -1; c = gc(); }

        T val = 0;
        while (c > ' ') {
            val = val * 10 + (c - '0');
            c = gc();
        }
        out = val * sign;
        return true;
    }
};

struct Node {
    int len;   // segment length (0 for padding)
    int pref;  // max prefix ones
    int suff;  // max suffix ones
    int best;  // max consecutive ones
};

static inline Node mergeNode(const Node &L, const Node &R) {
    if (L.len == 0) return R;
    if (R.len == 0) return L;
    Node t;
    t.len = L.len + R.len;
    t.pref = (L.pref == L.len) ? (L.len + R.pref) : L.pref;
    t.suff = (R.suff == R.len) ? (R.len + L.suff) : R.suff;
    t.best = max({L.best, R.best, L.suff + R.pref});
    return t;
}

struct SegTreeRuns {
    int n, size;
    vector<Node> seg;

    explicit SegTreeRuns(int n_) : n(n_) {
        size = 1;
        while (size < n) size <<= 1;
        seg.assign(2 * size, Node{0,0,0,0});

        // init lengths
        for (int i = 0; i < size; i++) {
            seg[size + i].len = (i < n) ? 1 : 0;
        }
        for (int i = size - 1; i >= 1; --i) {
            seg[i].len = seg[i<<1].len + seg[i<<1|1].len;
        }
    }

    void resetAll() {
        for (int i = 1; i < 2 * size; i++) {
            seg[i].pref = seg[i].suff = seg[i].best = 0;
        }
    }

    void activate(int idx) { // idx: 0-based, set to 1
        int p = size + idx;
        seg[p].pref = seg[p].suff = seg[p].best = 1;
        for (p >>= 1; p >= 1; p >>= 1) {
            seg[p] = mergeNode(seg[p<<1], seg[p<<1|1]);
        }
    }

    Node query(int l, int r) { // inclusive, 0-based
        Node left{0,0,0,0}, right{0,0,0,0};
        int L = l + size;
        int R = r + size + 1; // exclusive
        while (L < R) {
            if (L & 1) left = mergeNode(left, seg[L++]);
            if (R & 1) right = mergeNode(seg[--R], right);
            L >>= 1; R >>= 1;
        }
        return mergeNode(left, right);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    FastScanner fs;

    int N;
    fs.readInt(N);
    vector<long long> h(N);
    for (int i = 0; i < N; i++) fs.readInt(h[i]);

    // unique heights ascending
    vector<long long> vals = h;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    int Mv = (int)vals.size();

    // bars sorted by height desc
    vector<pair<long long,int>> bars;
    bars.reserve(N);
    for (int i = 0; i < N; i++) bars.push_back({h[i], i});
    sort(bars.begin(), bars.end(), [](auto &a, auto &b){
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });

    int Q;
    fs.readInt(Q);
    vector<int> L(Q), R(Q), W(Q);
    for (int i = 0; i < Q; i++) {
        int l, r, w;
        fs.readInt(l); fs.readInt(r); fs.readInt(w);
        --l; --r;
        L[i] = l; R[i] = r; W[i] = w;
    }

    vector<int> lo(Q, 0), hi(Q, Mv - 1);

    vector<vector<int>> bucket(Mv);
    SegTreeRuns seg(N);

    while (true) {
        bool any = false;
        for (int i = 0; i < Q; i++) {
            if (lo[i] < hi[i]) {
                int mid = (lo[i] + hi[i] + 1) >> 1;
                bucket[mid].push_back(i);
                any = true;
            }
        }
        if (!any) break;

        seg.resetAll();
        int p = 0; // pointer in bars

        // process mids from high to low (higher height -> fewer active)
        for (int mid = Mv - 1; mid >= 0; --mid) {
            if (bucket[mid].empty()) continue;

            long long thr = vals[mid];
            while (p < N && bars[p].first >= thr) {
                seg.activate(bars[p].second);
                p++;
            }

            for (int qi : bucket[mid]) {
                Node res = seg.query(L[qi], R[qi]);
                if (res.best >= W[qi]) lo[qi] = mid;
                else hi[qi] = mid - 1;
            }
            bucket[mid].clear();
        }
    }

    for (int i = 0; i < Q; i++) {
        cout << vals[lo[i]] << "\n";
    }
    return 0;
}
