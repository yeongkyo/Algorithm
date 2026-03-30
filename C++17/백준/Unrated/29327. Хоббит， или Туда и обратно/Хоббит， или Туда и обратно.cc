#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const long long INF = 1e18;

struct Point {
    long long x;
    int id;
    bool operator<(const Point& other) const {
        return x < other.x;
    }
};

int n;
vector<Point> pts;
vector<vector<long long>> C;

long long get_cost(int u, int v) {
    if (u == 0 || v == 0) return 0;
    return C[pts[u].id][pts[v].id];
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    pts.resize(n + 1);
    pts[0] = {0, 0};
    for (int i = 1; i <= n; ++i) {
        cin >> pts[i].x;
        pts[i].id = i;
    }

    sort(pts.begin() + 1, pts.end());

    C.assign(n + 1, vector<long long>(n + 1, 0));
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            cin >> C[i][j];
        }
    }

    vector<vector<long long>> dp(n + 1, vector<long long>(n + 1, INF));
    vector<vector<int>> prev_j(n + 1, vector<int>(n + 1, -1));

    dp[1][0] = 0;

    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (dp[i][j] == INF) continue;

            long long w1 = dp[i][j] + get_cost(i, i + 1);
            if (w1 < dp[i + 1][j]) {
                dp[i + 1][j] = w1;
                prev_j[i + 1][j] = j;
            }

            long long w2 = dp[i][j] + get_cost(j, i + 1);
            if (w2 < dp[i + 1][i]) {
                dp[i + 1][i] = w2;
                prev_j[i + 1][i] = j;
            }
        }
    }

    long long min_danger = INF;
    int opt_j = -1;

    for (int j = 0; j < n; ++j) {
        long long cur = dp[n][j] + get_cost(n, j);
        if (cur < min_danger) {
            min_danger = cur;
            opt_j = j;
        }
    }

    long long total_dist = 2LL * pts[n].x;

    vector<int> seq1 = {n};
    vector<int> seq2 = {opt_j};
    int curr_j = opt_j;

    for (int i = n; i >= 2; --i) {
        if (curr_j == i - 1) {
            int p = prev_j[i][i - 1];
            if (seq1.back() == i) seq1.push_back(p);
            else seq2.push_back(p);
            curr_j = p;
        } else {
            if (seq1.back() == i) seq1.push_back(i - 1);
            else seq2.push_back(i - 1);
        }
    }

    if (seq1.back() != 0) seq1.push_back(0);
    if (seq2.back() != 0) seq2.push_back(0);

    vector<int> final_path;
    for (int k = (int)seq2.size() - 1; k >= 0; --k) {
        if (seq2[k] != 0) final_path.push_back(pts[seq2[k]].id);
    }
    for (int k = 0; k < (int)seq1.size(); ++k) {
        if (seq1[k] != 0) final_path.push_back(pts[seq1[k]].id);
    }

    cout << total_dist << " " << min_danger << "\n";
    for (int i = 0; i < n; ++i) {
        cout << final_path[i] << (i == n - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}