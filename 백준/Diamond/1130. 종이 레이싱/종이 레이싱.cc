#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;

int N, M;
vector<string> grid;

bool visited[55][55][131][131];

int check_path(int r1, int c1, int r2, int c2) {
    int dr = r2 - r1;
    int dc = c2 - c1;
    if (dr == 0 && dc == 0) {
        if (r1 < 0 || r1 >= N || c1 < 0 || c1 >= M) return 0;
        if (grid[r1][c1] == 'F') return 2;
        if (grid[r1][c1] == 'X') return 0;
        return 1;
    }

    vector<double> ts;
    ts.push_back(1.0);

    if (dr != 0) {
        int min_r = min(r1, r2);
        int max_r = max(r1, r2);
        for (int x = min_r; x < max_r; ++x) {
            double t = (x + 0.5 - r1) / (double)dr;
            if (t > 1e-9 && t < 1.0 - 1e-9) ts.push_back(t);
        }
    }
    if (dc != 0) {
        int min_c = min(c1, c2);
        int max_c = max(c1, c2);
        for (int y = min_c; y < max_c; ++y) {
            double t = (y + 0.5 - c1) / (double)dc;
            if (t > 1e-9 && t < 1.0 - 1e-9) ts.push_back(t);
        }
    }

    sort(ts.begin(), ts.end());
    vector<double> unique_ts;
    for (double t : ts) {
        if (unique_ts.empty() || t - unique_ts.back() > 1e-9) {
            unique_ts.push_back(t);
        }
    }

    double prev_t = 0.0;
    for (double t : unique_ts) {
        double t_mid = (prev_t + t) / 2.0;
        double r_mid = r1 + t_mid * dr;
        double c_mid = c1 + t_mid * dc;
        int cr = round(r_mid);
        int cc = round(c_mid);

        if (cr < 0 || cr >= N || cc < 0 || cc >= M) return 0;
        if (grid[cr][cc] == 'X') return 0;
        if (grid[cr][cc] == 'F') return 2;

        double r_pt = r1 + t * dr;
        double c_pt = c1 + t * dc;

        int min_cr = ceil(r_pt - 0.5 - 1e-7);
        int max_cr = floor(r_pt + 0.5 + 1e-7);
        int min_cc = ceil(c_pt - 0.5 - 1e-7);
        int max_cc = floor(c_pt + 0.5 + 1e-7);

        for (int i = min_cr; i <= max_cr; ++i) {
            for (int j = min_cc; j <= max_cc; ++j) {
                if (i >= 0 && i < N && j >= 0 && j < M) {
                    if (grid[i][j] == 'F') return 2;
                }
            }
        }

        prev_t = t;
    }

    return 1;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M)) return 0;

    grid.resize(N);
    int sr = -1, sc = -1;
    for (int i = 0; i < N; ++i) {
        cin >> grid[i];
        for (int j = 0; j < M; ++j) {
            if (grid[i][j] == 'S') {
                sr = i;
                sc = j;
            }
        }
    }

    int vr_init, vc_init;
    cin >> vr_init >> vc_init;

    queue<tuple<int, int, int, int, int>> q;
    q.push({sr, sc, vr_init, vc_init, 0});
    visited[sr][sc][vr_init + 65][vc_init + 65] = true;

    while (!q.empty()) {
        auto [r, c, vr, vc, turns] = q.front();
        q.pop();

        for (int dvr = -1; dvr <= 1; ++dvr) {
            for (int dvc = -1; dvc <= 1; ++dvc) {
                int nvr = vr + dvr;
                int nvc = vc + dvc;
                int nr = r + nvr;
                int nc = c + nvc;

                int res = check_path(r, c, nr, nc);
                if (res == 2) {
                    cout << turns + 1 << "\n";
                    return 0;
                }
                if (res == 1) {
                    if (nvr + 65 >= 0 && nvr + 65 < 131 && nvc + 65 >= 0 && nvc + 65 < 131) {
                        if (!visited[nr][nc][nvr + 65][nvc + 65]) {
                            visited[nr][nc][nvr + 65][nvc + 65] = true;
                            q.push({nr, nc, nvr, nvc, turns + 1});
                        }
                    }
                }
            }
        }
    }

    cout << -1 << "\n";
    return 0;
}