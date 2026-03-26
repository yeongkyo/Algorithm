#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

struct Interval {
    int S, E, id;
    long long V;
    bool operator<(const Interval& o) const {
        if (S != o.S) return S < o.S;
        return E > o.E;
    }
};

mt19937 rng(1337);

struct Node {
    long long val;
    int pri, sz;
    Node *l, *r;
    Node(long long v) : val(v), pri(rng()), sz(1), l(nullptr), r(nullptr) {}
};

int sz(Node* t) { return t ? t->sz : 0; }
void upd(Node* t) { if (t) t->sz = 1 + sz(t->l) + sz(t->r); }

void split(Node* t, int k, Node*& l, Node*& r) {
    if (!t) { l = r = nullptr; return; }
    if (sz(t->l) >= k) {
        split(t->l, k, l, t->l);
        r = t;
    } else {
        split(t->r, k - sz(t->l) - 1, t->r, r);
        l = t;
    }
    upd(t);
}

void merge(Node*& t, Node* l, Node* r) {
    if (!l || !r) { t = l ? l : r; return; }
    if (l->pri > r->pri) {
        merge(l->r, l->r, r);
        t = l;
    } else {
        merge(r->l, l, r->l);
        t = r;
    }
    upd(t);
}

int get_pos(Node* t, long long v) {
    if (!t) return 0;
    if (t->val < v) {
        return get_pos(t->l, v);
    } else {
        return sz(t->l) + 1 + get_pos(t->r, v);
    }
}

void extract_all(Node* t, vector<long long>& A) {
    if (!t) return;
    extract_all(t->l, A);
    A.push_back(t->val);
    extract_all(t->r, A);
}

void add_to(Node* t, const vector<long long>& A, int& idx) {
    if (!t) return;
    add_to(t->l, A, idx);
    t->val += A[idx++];
    add_to(t->r, A, idx);
}

int N;
vector<Interval> intervals;
vector<vector<int>> adj;

Node* dfs(int u) {
    Node* L = nullptr;
    for(int v : adj[u]) {
        Node* S = dfs(v);
        if (!L) {
            L = S;
        } else {
            if (sz(L) < sz(S)) swap(L, S);
            if (sz(S) > 0) {
                vector<long long> A;
                extract_all(S, A);
                Node *L1, *L2;
                split(L, A.size(), L1, L2);
                int idx = 0;
                add_to(L1, A, idx);
                merge(L, L1, L2);
            }
        }
    }
    long long V = intervals[u].V;
    int pos = get_pos(L, V);
    Node *L1, *L2;
    split(L, pos, L1, L2);
    Node *mid = new Node(V);
    merge(L, L1, mid);
    merge(L, L, L2);
    return L;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    if (!(cin >> N)) return 0;
    intervals.resize(N);
    for(int i=0; i<N; ++i) {
        cin >> intervals[i].S >> intervals[i].E >> intervals[i].V;
        intervals[i].id = i;
    }
    sort(intervals.begin(), intervals.end());

    vector<int> roots;
    adj.resize(N);
    vector<int> st;

    for(int i=0; i<N; ++i) {
        while(!st.empty() && intervals[st.back()].E <= intervals[i].S) {
            st.pop_back();
        }
        if(st.empty()) {
            roots.push_back(i);
        } else {
            adj[st.back()].push_back(i);
        }
        st.push_back(i);
    }

    Node* global_L = nullptr;
    for(int r : roots) {
        Node* S = dfs(r);
        if (!global_L) {
            global_L = S;
        } else {
            if (sz(global_L) < sz(S)) swap(global_L, S);
            if (sz(S) > 0) {
                vector<long long> A;
                extract_all(S, A);
                Node *L1, *L2;
                split(global_L, A.size(), L1, L2);
                int idx = 0;
                add_to(L1, A, idx);
                merge(global_L, L1, L2);
            }
        }
    }

    vector<long long> final_D;
    extract_all(global_L, final_D);

    long long sum = 0;
    for(int k=1; k<=N; ++k) {
        if(k - 1 < (int)final_D.size()) {
            sum += final_D[k - 1];
        }
        cout << sum << (k == N ? "" : " ");
    }
    cout << "\n";
    return 0;
}