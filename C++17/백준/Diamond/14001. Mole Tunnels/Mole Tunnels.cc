#include <iostream>
#include <vector>

using namespace std;

const long long INF = 1e15;

int n, m;
vector<int> c;
vector<int> p_arr;
vector<int> flow_arr;
vector<pair<long long, int>> min_dist;

void update(int u) {
    while (u > 0) {
        pair<long long, int> res = {INF, -1};
        if (c[u] > 0) res = {0, u};
        
        if (2 * u <= n) {
            int w = 2 * u;
            if (min_dist[w].first != INF) {
                long long cost = (flow_arr[w] > 0 ? -1 : 1);
                if (cost + min_dist[w].first < res.first) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                } else if (cost + min_dist[w].first == res.first && min_dist[w].second < res.second) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                }
            }
        }
        if (2 * u + 1 <= n) {
            int w = 2 * u + 1;
            if (min_dist[w].first != INF) {
                long long cost = (flow_arr[w] > 0 ? -1 : 1);
                if (cost + min_dist[w].first < res.first) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                } else if (cost + min_dist[w].first == res.first && min_dist[w].second < res.second) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                }
            }
        }
        min_dist[u] = res;
        u /= 2;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n >> m)) return 0;

    c.resize(n + 1);
    for (int i = 1; i <= n; ++i) {
        cin >> c[i];
    }

    p_arr.resize(m + 1);
    for (int i = 1; i <= m; ++i) {
        cin >> p_arr[i];
    }

    flow_arr.assign(n + 1, 0);
    min_dist.assign(n + 1, {INF, -1});

    for (int u = n; u >= 1; --u) {
        pair<long long, int> res = {INF, -1};
        if (c[u] > 0) res = {0, u};
        
        if (2 * u <= n) {
            int w = 2 * u;
            if (min_dist[w].first != INF) {
                long long cost = 1; 
                if (cost + min_dist[w].first < res.first) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                } else if (cost + min_dist[w].first == res.first && min_dist[w].second < res.second) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                }
            }
        }
        if (2 * u + 1 <= n) {
            int w = 2 * u + 1;
            if (min_dist[w].first != INF) {
                long long cost = 1;
                if (cost + min_dist[w].first < res.first) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                } else if (cost + min_dist[w].first == res.first && min_dist[w].second < res.second) {
                    res = {cost + min_dist[w].first, min_dist[w].second};
                }
            }
        }
        min_dist[u] = res;
    }

    long long total_ans = 0;
    for (int k = 1; k <= m; ++k) {
        int p = p_arr[k];
        
        long long best_total = INF;
        int best_A = -1;
        int best_v = -1;
        
        int curr_A = p;
        long long up_cost = 0;
        
        while (curr_A > 0) {
            if (min_dist[curr_A].first != INF) {
                if (up_cost + min_dist[curr_A].first < best_total) {
                    best_total = up_cost + min_dist[curr_A].first;
                    best_A = curr_A;
                    best_v = min_dist[curr_A].second;
                } else if (up_cost + min_dist[curr_A].first == best_total) {
                    if (min_dist[curr_A].second < best_v) {
                        best_A = curr_A;
                        best_v = min_dist[curr_A].second;
                    }
                }
            }
            if (curr_A == 1) break;
            long long cost_up = (flow_arr[curr_A] < 0 ? -1 : 1);
            up_cost += cost_up;
            curr_A /= 2;
        }
        
        total_ans += best_total;
        cout << total_ans << (k == m ? "" : " ");
        
        c[best_v]--;
        
        int curr = p;
        while (curr != best_A) {
            flow_arr[curr]++;
            curr /= 2;
        }
        
        curr = best_v;
        while (curr != best_A) {
            flow_arr[curr]--;
            curr /= 2;
        }
        
        update(p);
        update(best_v);
    }
    cout << "\n";

    return 0;
}