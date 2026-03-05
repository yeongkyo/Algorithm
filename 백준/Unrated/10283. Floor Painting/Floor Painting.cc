#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct ExtRect {
    int left, right, bottom, top;
};

struct Edge {
    int x, y1, y2;
};

struct Event {
    int x, type, y1, y2;
};

struct ActiveRect {
    int left, right, bottom;
};

const int MAXM = 100005;
int tree_arr[4 * MAXM];
int lazy_arr[4 * MAXM];

void push(int node) {
    if (lazy_arr[node]) {
        tree_arr[2 * node] += lazy_arr[node];
        lazy_arr[2 * node] += lazy_arr[node];
        tree_arr[2 * node + 1] += lazy_arr[node];
        lazy_arr[2 * node + 1] += lazy_arr[node];
        lazy_arr[node] = 0;
    }
}

void update(int node, int l, int r, int ql, int qr, int val) {
    if (ql <= l && r <= qr) {
        tree_arr[node] += val;
        lazy_arr[node] += val;
        return;
    }
    push(node);
    int mid = (l + r) / 2;
    if (ql <= mid) update(2 * node, l, mid, ql, qr, val);
    if (qr > mid) update(2 * node + 1, mid + 1, r, ql, qr, val);
    tree_arr[node] = min(tree_arr[2 * node], tree_arr[2 * node + 1]);
}

bool check(int S, const vector<ExtRect>& ext_rects, int X_min, int X_max, int Y_min, int Y_max) {
    int M = Y_max - S - Y_min;
    if (M < 0) return false;
    if (X_max - S < X_min) return false;

    vector<Event> events;
    for (const auto& r : ext_rects) {
        int L = max(X_min, r.left - S + 1);
        int R = r.right;
        if (L >= R) continue;

        int bottom = max(Y_min, r.bottom - S + 1) - Y_min;
        int top = min(Y_max - S, r.top - 1) - Y_min;
        if (bottom > top) continue;

        events.push_back({L, 1, bottom, top});
        events.push_back({R, -1, bottom, top});
    }
    events.push_back({X_min, 0, 0, 0}); 

    sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        return a.x < b.x;
    });

    for (int i = 0; i < 4 * (M + 2); ++i) {
        tree_arr[i] = 0;
        lazy_arr[i] = 0;
    }

    int i = 0;
    while (i < events.size()) {
        int curr_x = events[i].x;
        if (curr_x > X_max - S) break;

        while (i < events.size() && events[i].x == curr_x) {
            if (events[i].type != 0) {
                update(1, 0, M, events[i].y1, events[i].y2, events[i].type);
            }
            i++;
        }

        if (curr_x >= X_min && tree_arr[1] == 0) {
            return true;
        }
    }
    return false;
}

void solve() {
    int n;
    if (!(cin >> n)) return;

    vector<int> X(n), Y(n);
    int X_min = 1e9, X_max = -1e9;
    int Y_min = 1e9, Y_max = -1e9;

    for (int i = 0; i < n; ++i) {
        cin >> X[i] >> Y[i];
        X_min = min(X_min, X[i]);
        X_max = max(X_max, X[i]);
        Y_min = min(Y_min, Y[i]);
        Y_max = max(Y_max, Y[i]);
    }

    int BX_min = X_min - 1, BX_max = X_max + 1;
    int BY_min = Y_min - 1, BY_max = Y_max + 1;

    vector<Edge> vertical_edges;
    for (int i = 0; i < n; ++i) {
        int nx = (i == n - 1) ? 0 : i + 1;
        if (X[i] == X[nx]) {
            vertical_edges.push_back({X[i], min(Y[i], Y[nx]), max(Y[i], Y[nx])});
        }
    }

    vector<int> unique_y = Y;
    sort(unique_y.begin(), unique_y.end());
    unique_y.erase(unique(unique_y.begin(), unique_y.end()), unique_y.end());

    vector<ExtRect> ext_rects;
    vector<ActiveRect> current_rects;

    for (size_t j = 0; j <= unique_y.size(); ++j) {
        int b_bottom = (j == 0) ? BY_min : unique_y[j - 1];
        int b_top = (j == unique_y.size()) ? BY_max : unique_y[j];
        double y_mid = b_bottom + (b_top - b_bottom) / 2.0;

        vector<int> active_x;
        for (const auto& edge : vertical_edges) {
            if (edge.y1 < y_mid && edge.y2 > y_mid) {
                active_x.push_back(edge.x);
            }
        }
        sort(active_x.begin(), active_x.end());

        vector<pair<int, int>> new_intervals;
        int curr_x = BX_min;
        for (size_t i = 0; i < active_x.size(); i += 2) {
            if (curr_x < active_x[i]) {
                new_intervals.push_back({curr_x, active_x[i]});
            }
            curr_x = active_x[i + 1];
        }
        if (curr_x < BX_max) {
            new_intervals.push_back({curr_x, BX_max});
        }

        vector<ActiveRect> next_rects;
        for (const auto& n_intv : new_intervals) {
            bool matched = false;
            for (const auto& c_rect : current_rects) {
                if (n_intv.first == c_rect.left && n_intv.second == c_rect.right) {
                    next_rects.push_back({c_rect.left, c_rect.right, c_rect.bottom});
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                next_rects.push_back({n_intv.first, n_intv.second, b_bottom});
            }
        }

        for (const auto& c_rect : current_rects) {
            bool matched = false;
            for (const auto& n_intv : new_intervals) {
                if (n_intv.first == c_rect.left && n_intv.second == c_rect.right) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                ext_rects.push_back({c_rect.left, c_rect.right, c_rect.bottom, b_bottom});
            }
        }
        current_rects = next_rects;
    }

    for (const auto& c_rect : current_rects) {
        ext_rects.push_back({c_rect.left, c_rect.right, c_rect.bottom, BY_max});
    }

    int ans = 0;
    int low = 1, high = min(X_max - X_min, Y_max - Y_min);
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (check(mid, ext_rects, X_min, X_max, Y_min, Y_max)) {
            ans = mid;
            low = mid + 1; 
        } else {
            high = mid - 1;
        }
    }

    cout << ans << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t;
    if (cin >> t) {
        while (t--) {
            solve();
        }
    }
    return 0;
}