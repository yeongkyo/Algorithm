#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int v_x[20005], v_ymin[20005], v_ymax[20005];
int h_y[20005], h_xmin[20005], h_xmax[20005];
int v_cnt = 0, h_cnt = 0;

struct Point { int x, y; };
vector<Point> targets;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> vx(n), vy(n);
    for (int i = 0; i < n; ++i) {
        cin >> vx[i] >> vy[i];
        vx[i] *= 4;
        vy[i] *= 4;
    }

    for (int i = 0; i < n; ++i) {
        int nx = vx[(i + 1) % n];
        int ny = vy[(i + 1) % n];
        
        if (vx[i] == nx) {
            v_x[v_cnt] = vx[i];
            v_ymin[v_cnt] = min(vy[i], ny);
            v_ymax[v_cnt] = max(vy[i], ny);
            v_cnt++;
        } else {
            h_y[h_cnt] = vy[i];
            h_xmin[h_cnt] = min(vx[i], nx);
            h_xmax[h_cnt] = max(vx[i], nx);
            h_cnt++;
        }

        int prev = (i - 1 + n) % n;
        long long dx_in = vx[i] - vx[prev];
        long long dy_in = vy[i] - vy[prev];
        long long dx_out = nx - vx[i];
        long long dy_out = ny - vy[i];
        
        int dxi = (dx_in > 0) ? 1 : ((dx_in < 0) ? -1 : 0);
        int dyi = (dy_in > 0) ? 1 : ((dy_in < 0) ? -1 : 0);
        int dxo = (dx_out > 0) ? 1 : ((dx_out < 0) ? -1 : 0);
        int dyo = (dy_out > 0) ? 1 : ((dy_out < 0) ? -1 : 0);
        
        if (dxi * dyo - dyi * dxo > 0) {
            targets.push_back({vx[i] - dyi - dyo, vy[i] + dxi + dxo});
        }
    }

    for (int q = 0; q < 5; ++q) {
        int qx, qy;
        cin >> qx >> qy;
        int fx = qx * 4;
        int fy = qy * 4;

        bool all_reachable = true;

        for (const auto& P : targets) {
            bool p1_ok = true;
            int min_x = min(fx, P.x), max_x = max(fx, P.x);
            for (int i = 0; i < v_cnt; ++i) {
                if (v_x[i] > min_x && v_x[i] < max_x && v_ymin[i] < fy && fy < v_ymax[i]) {
                    p1_ok = false; 
                    break;
                }
            }
            if (p1_ok) {
                int min_y = min(fy, P.y), max_y = max(fy, P.y);
                for (int i = 0; i < h_cnt; ++i) {
                    if (h_y[i] > min_y && h_y[i] < max_y && h_xmin[i] < P.x && P.x < h_xmax[i]) {
                        p1_ok = false; 
                        break;
                    }
                }
            }

            if (p1_ok) continue;

            bool p2_ok = true;
            int min_y2 = min(fy, P.y), max_y2 = max(fy, P.y);
            for (int i = 0; i < h_cnt; ++i) {
                if (h_y[i] > min_y2 && h_y[i] < max_y2 && h_xmin[i] < fx && fx < h_xmax[i]) {
                    p2_ok = false; 
                    break;
                }
            }
            if (p2_ok) {
                int min_x2 = min(fx, P.x), max_x2 = max(fx, P.x);
                for (int i = 0; i < v_cnt; ++i) {
                    if (v_x[i] > min_x2 && v_x[i] < max_x2 && v_ymin[i] < P.y && P.y < v_ymax[i]) {
                        p2_ok = false; 
                        break;
                    }
                }
            }

            if (!p2_ok) {
                all_reachable = false;
                break;
            }
        }

        if (all_reachable) {
            cout << "YES\n";
        } else {
            cout << "NO\n";
        }
    }

    return 0;
}