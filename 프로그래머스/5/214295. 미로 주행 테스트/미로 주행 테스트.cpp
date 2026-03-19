#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MAX_S = 1000005;
const long long INF = 2e18;

struct Node {
    int cov;
    long long len0, len1;
} tree[4 * MAX_S];

long long pre_len0[MAX_S];
long long pre_len1[MAX_S];
vector<long long> y_vals;

inline long long count_parity(long long L, long long R, int parity) {
    if (L > R) return 0;
    long long first_val = L;
    if ((first_val & 1) != parity) first_val++;
    if (first_val > R) return 0;
    long long last_val = R;
    if ((last_val & 1) != parity) last_val--;
    return (last_val - first_val) / 2 + 1;
}

void push_up(int node, int l, int r) {
    if (tree[node].cov > 0) {
        tree[node].len0 = 0;
        tree[node].len1 = 0;
    } else {
        if (l == r) {
            tree[node].len0 = pre_len0[l];
            tree[node].len1 = pre_len1[l];
        } else {
            tree[node].len0 = tree[2*node].len0 + tree[2*node+1].len0;
            tree[node].len1 = tree[2*node].len1 + tree[2*node+1].len1;
        }
    }
}

void build(int node, int l, int r) {
    tree[node].cov = 0;
    if (l == r) {
        tree[node].len0 = pre_len0[l];
        tree[node].len1 = pre_len1[l];
        return;
    }
    int mid = l + (r - l) / 2;
    build(2*node, l, mid);
    build(2*node+1, mid+1, r);
    push_up(node, l, r);
}

void update(int node, int l, int r, int ql, int qr, int val) {
    if (ql > r || qr < l) return;
    if (ql <= l && r <= qr) {
        tree[node].cov += val;
        push_up(node, l, r);
        return;
    }
    int mid = l + (r - l) / 2;
    update(2*node, l, mid, ql, qr, val);
    update(2*node+1, mid+1, r, ql, qr, val);
    push_up(node, l, r);
}

pair<long long, long long> query(int node, int l, int r, long long ql_val, long long qr_val, int cov_sum) {
    long long node_l_val = y_vals[l];
    long long node_r_val = y_vals[r+1] - 1;
    
    if (ql_val > node_r_val || qr_val < node_l_val) return {0, 0};
    
    if (cov_sum > 0 || tree[node].cov > 0) {
        return {0, 0};
    }
    
    if (ql_val <= node_l_val && node_r_val <= qr_val) {
        return {tree[node].len0, tree[node].len1};
    }
    
    if (l == r) {
        long long act_l = max(node_l_val, ql_val);
        long long act_r = min(node_r_val, qr_val);
        if (act_l > act_r) return {0, 0};
        long long len0 = count_parity(act_l, act_r, 0);
        long long len1 = count_parity(act_l, act_r, 1);
        return {len0, len1};
    }
    
    int mid = l + (r - l) / 2;
    auto left_res = query(2*node, l, mid, ql_val, qr_val, cov_sum + tree[node].cov);
    auto right_res = query(2*node+1, mid+1, r, ql_val, qr_val, cov_sum + tree[node].cov);
    return {left_res.first + right_res.first, left_res.second + right_res.second};
}

struct Rect {
    long long u_min, u_max;
    long long v_min, v_max;
};

struct Event {
    long long u;
    int type;
    long long v_min, v_max_plus_1;
};

long long solution(int n, int m, vector<vector<int>> tests) {
    long long U_min = 0;
    long long U_max = (long long)n + m;
    long long V_min = -(long long)m;
    long long V_max = n;

    long long U1_min = -INF, U1_max = INF;
    long long V1_min = -INF, V1_max = INF;

    vector<Rect> rects;
    rects.reserve(tests.size());

    for (const auto& t : tests) {
        long long x = t[0], y = t[1], d = t[2], flag = t[3];
        long long u0 = x + y, v0 = x - y;
        if (flag == 1) {
            U1_min = max(U1_min, u0 - d);
            U1_max = min(U1_max, u0 + d);
            V1_min = max(V1_min, v0 - d);
            V1_max = min(V1_max, v0 + d);
        } else {
            rects.push_back({u0 - d, u0 + d, v0 - d, v0 + d});
        }
    }

    U_min = max(U_min, U1_min);
    U_max = min(U_max, U1_max);
    V_min = max(V_min, V1_min);
    V_max = min(V_max, V1_max);

    if (U_min > U_max || V_min > V_max) return 0;

    y_vals.clear();
    y_vals.push_back(V_min);
    y_vals.push_back(V_max + 1);
    for (const auto& r : rects) {
        y_vals.push_back(r.v_min);
        y_vals.push_back(r.v_max + 1);
    }
    sort(y_vals.begin(), y_vals.end());
    y_vals.erase(unique(y_vals.begin(), y_vals.end()), y_vals.end());

    int S = y_vals.size() - 1;
    for (int i = 0; i < S; ++i) {
        pre_len0[i] = count_parity(y_vals[i], y_vals[i+1] - 1, 0);
        pre_len1[i] = count_parity(y_vals[i], y_vals[i+1] - 1, 1);
    }
    
    build(1, 0, S - 1);

    vector<long long> U_events;
    U_events.reserve(4000000); 
    auto add_U = [&](long long u) {
        if (u >= U_min && u <= U_max + 1) {
            U_events.push_back(u);
        }
    };

    auto add_V_events = [&](long long v) {
        add_U(-v);
        add_U(v + 2LL * m);
        add_U(v);
        add_U(2LL * n - v);
    };

    add_U(U_min); add_U(U_max + 1);
    for (const auto& r : rects) {
        add_U(r.u_min); add_U(r.u_max + 1);
    }

    add_U(m); add_U(n); add_U(0); add_U((long long)n + m);

    for (long long v : y_vals) {
        add_V_events(v);
        add_V_events(v - 1);
    }

    sort(U_events.begin(), U_events.end());
    U_events.erase(unique(U_events.begin(), U_events.end()), U_events.end());

    vector<Event> events;
    events.reserve(rects.size() * 2);
    for (const auto& r : rects) {
        events.push_back({r.u_min, 1, r.v_min, r.v_max + 1});
        events.push_back({r.u_max + 1, -1, r.v_min, r.v_max + 1});
    }
    sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        return a.u < b.u;
    });

    int event_idx = 0;
    long long total_valid = 0;

    auto get_g = [&](long long u) -> long long {
        long long l_val = max({V_min, -u, u - 2LL * m});
        long long r_val = min({V_max, u, 2LL * n - u});
        if (l_val > r_val) return 0;
        
        auto res = query(1, 0, S - 1, l_val, r_val, 0);
        return ((u & 1) == 0) ? res.first : res.second;
    };

    for (size_t k = 0; k + 1 < U_events.size(); ++k) {
        long long curr_u = U_events[k];
        long long next_u = U_events[k+1];
        
        while (event_idx < events.size() && events[event_idx].u <= curr_u) {
            int ql = lower_bound(y_vals.begin(), y_vals.end(), events[event_idx].v_min) - y_vals.begin();
            int qr = lower_bound(y_vals.begin(), y_vals.end(), events[event_idx].v_max_plus_1) - y_vals.begin() - 1;
            if (ql <= qr) {
                update(1, 0, S - 1, ql, qr, events[event_idx].type);
            }
            event_idx++;
        }
        
        long long len = next_u - curr_u;
        if (len == 0) continue;
        
        if (len <= 8) {
            for (long long u = curr_u; u < next_u; ++u) total_valid += get_g(u);
        } else {
            long long u0 = curr_u, u1 = curr_u + 1, u2 = curr_u + 2, u3 = curr_u + 3;
            
            long long g0 = get_g(u0);
            long long g2 = get_g(u2);
            long long diff0 = g2 - g0;
            
            long long g1 = get_g(u1);
            long long g3 = get_g(u3);
            long long diff1 = g3 - g1;
            
            long long N0 = (next_u - 1 - u0) / 2 + 1;
            long long N1 = (next_u - 1 - u1) / 2 + 1;
            
            long long pairs0 = (long long)N0 * (N0 - 1) / 2;
            long long sum0 = (long long)N0 * g0 + pairs0 * diff0;
            
            long long pairs1 = (long long)N1 * (N1 - 1) / 2;
            long long sum1 = (long long)N1 * g1 + pairs1 * diff1;
            
            total_valid += sum0 + sum1;
        }
    }

    return total_valid;
}