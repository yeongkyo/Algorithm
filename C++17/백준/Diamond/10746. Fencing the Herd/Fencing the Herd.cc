#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const long long INF = 4000000000000000000LL;

struct Point {
    long long x, y;
    bool operator<(const Point& o) const {
        if (x != o.x) return x < o.x;
        return y < o.y;
    }
};

__int128 cross(Point p1, Point p2, Point p3) {
    return (__int128)(p2.x - p1.x) * (p3.y - p1.y) - (__int128)(p2.y - p1.y) * (p3.x - p1.x);
}

long long eval(Point p, long long A, long long B) {
    return A * p.x + B * p.y;
}

long long get_max_upper(const vector<Point>& hull, long long A, long long B) {
    if (hull.empty()) return -INF;
    if (B >= 0) {
        int l = 0, r = hull.size() - 1;
        while (r - l >= 3) {
            int m1 = l + (r - l) / 3;
            int m2 = r - (r - l) / 3;
            if (eval(hull[m1], A, B) < eval(hull[m2], A, B)) l = m1;
            else r = m2;
        }
        long long mx = -INF;
        for (int i = l; i <= r; ++i) mx = max(mx, eval(hull[i], A, B));
        return mx;
    } else {
        return max(eval(hull.front(), A, B), eval(hull.back(), A, B));
    }
}

long long get_min_upper(const vector<Point>& hull, long long A, long long B) {
    if (hull.empty()) return INF;
    if (B <= 0) {
        int l = 0, r = hull.size() - 1;
        while (r - l >= 3) {
            int m1 = l + (r - l) / 3;
            int m2 = r - (r - l) / 3;
            if (eval(hull[m1], A, B) > eval(hull[m2], A, B)) l = m1;
            else r = m2;
        }
        long long mn = INF;
        for (int i = l; i <= r; ++i) mn = min(mn, eval(hull[i], A, B));
        return mn;
    } else {
        return min(eval(hull.front(), A, B), eval(hull.back(), A, B));
    }
}

long long get_max_lower(const vector<Point>& hull, long long A, long long B) {
    if (hull.empty()) return -INF;
    if (B <= 0) {
        int l = 0, r = hull.size() - 1;
        while (r - l >= 3) {
            int m1 = l + (r - l) / 3;
            int m2 = r - (r - l) / 3;
            if (eval(hull[m1], A, B) < eval(hull[m2], A, B)) l = m1;
            else r = m2;
        }
        long long mx = -INF;
        for (int i = l; i <= r; ++i) mx = max(mx, eval(hull[i], A, B));
        return mx;
    } else {
        return max(eval(hull.front(), A, B), eval(hull.back(), A, B));
    }
}

long long get_min_lower(const vector<Point>& hull, long long A, long long B) {
    if (hull.empty()) return INF;
    if (B >= 0) {
        int l = 0, r = hull.size() - 1;
        while (r - l >= 3) {
            int m1 = l + (r - l) / 3;
            int m2 = r - (r - l) / 3;
            if (eval(hull[m1], A, B) > eval(hull[m2], A, B)) l = m1;
            else r = m2;
        }
        long long mn = INF;
        for (int i = l; i <= r; ++i) mn = min(mn, eval(hull[i], A, B));
        return mn;
    } else {
        return min(eval(hull.front(), A, B), eval(hull.back(), A, B));
    }
}

struct Node {
    vector<Point> upper, lower;
};

const int MAX_PTS = 200005;
Node tree[4 * MAX_PTS];
vector<Point> pts;

vector<Point> build(int node, int l, int r) {
    if (l == r) {
        tree[node].upper = {pts[l]};
        tree[node].lower = {pts[l]};
        return {pts[l]};
    }
    int mid = l + (r - l) / 2;
    vector<Point> left_pts = build(2 * node, l, mid);
    vector<Point> right_pts = build(2 * node + 1, mid + 1, r);
    
    vector<Point> sorted_pts;
    sorted_pts.reserve(left_pts.size() + right_pts.size());
    merge(left_pts.begin(), left_pts.end(),
          right_pts.begin(), right_pts.end(),
          back_inserter(sorted_pts));
    
    for (const auto& p : sorted_pts) {
        while (tree[node].upper.size() >= 2 && cross(tree[node].upper[tree[node].upper.size() - 2], tree[node].upper.back(), p) >= 0) {
            tree[node].upper.pop_back();
        }
        tree[node].upper.push_back(p);
    }
    
    for (const auto& p : sorted_pts) {
        while (tree[node].lower.size() >= 2 && cross(tree[node].lower[tree[node].lower.size() - 2], tree[node].lower.back(), p) <= 0) {
            tree[node].lower.pop_back();
        }
        tree[node].lower.push_back(p);
    }
    
    tree[node].upper.shrink_to_fit();
    tree[node].lower.shrink_to_fit();
    
    return sorted_pts;
}

long long global_max, global_min;

void query_tree(int node, int l, int r, int ql, int qr, long long A, long long B) {
    if (ql <= l && r <= qr) {
        global_max = max(global_max, get_max_upper(tree[node].upper, A, B));
        global_max = max(global_max, get_max_lower(tree[node].lower, A, B));
        global_min = min(global_min, get_min_upper(tree[node].upper, A, B));
        global_min = min(global_min, get_min_lower(tree[node].lower, A, B));
        return;
    }
    int mid = l + (r - l) / 2;
    if (ql <= mid) query_tree(2 * node, l, mid, ql, qr, A, B);
    if (qr > mid) query_tree(2 * node + 1, mid + 1, r, ql, qr, A, B);
}

struct Query {
    int type;
    long long A, B, C;
    int pt_idx;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n, q;
    if (!(cin >> n >> q)) return 0;
    
    pts.resize(n);
    for (int i = 0; i < n; ++i) {
        cin >> pts[i].x >> pts[i].y;
    }
    
    vector<Query> queries;
    queries.reserve(q);
    
    for (int i = 0; i < q; ++i) {
        int type;
        cin >> type;
        if (type == 1) {
            long long x, y;
            cin >> x >> y;
            pts.push_back({x, y});
        } else {
            long long A, B, C;
            cin >> A >> B >> C;
            queries.push_back({2, A, B, C, (int)pts.size() - 1});
        }
    }
    
    int S = pts.size();
    build(1, 0, S - 1);
    
    for (const auto& qry : queries) {
        global_max = -INF;
        global_min = INF;
        query_tree(1, 0, S - 1, 0, qry.pt_idx, qry.A, qry.B);
        
        if (global_max < qry.C || global_min > qry.C) {
            cout << "YES\n";
        } else {
            cout << "NO\n";
        }
    }
    
    return 0;
}