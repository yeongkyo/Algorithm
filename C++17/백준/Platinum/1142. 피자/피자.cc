#include <iostream>
#include <vector>

using namespace std;

struct Point {
    long long x, y;
};

long long dist2(Point p) {
    return p.x * p.x + p.y * p.y;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<Point> P(N);
    for (int i = 0; i < N; ++i) {
        cin >> P[i].x >> P[i].y;
    }

    int p0_idx = -1;
    for (int i = 0; i < N; ++i) {
        if (P[i].x != 0 || P[i].y != 0) {
            p0_idx = i;
            break;
        }
    }

    if (p0_idx == -1) {
        cout << -1 << "\n";
        return 0;
    }

    int ans = 0;
    for (int j = 0; j < N; ++j) {
        if (dist2(P[p0_idx]) != dist2(P[j])) continue;

        long long vx, vy;
        if (P[p0_idx].x == P[j].x && P[p0_idx].y == P[j].y) {
            vx = P[p0_idx].x;
            vy = P[p0_idx].y;
        } else if (P[p0_idx].x == -P[j].x && P[p0_idx].y == -P[j].y) {
            vx = -P[p0_idx].y;
            vy = P[p0_idx].x;
        } else {
            vx = P[p0_idx].x + P[j].x;
            vy = P[p0_idx].y + P[j].y;
        }

        long long D = vx * vx + vy * vy;
        if (D == 0) continue;

        bool valid = true;
        for (int k = 0; k < N; ++k) {
            long long dot = P[k].x * vx + P[k].y * vy;
            long long target_x = 2 * dot * vx - P[k].x * D;
            long long target_y = 2 * dot * vy - P[k].y * D;

            bool found = false;
            for (int m = 0; m < N; ++m) {
                if (P[m].x * D == target_x && P[m].y * D == target_y) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                valid = false;
                break;
            }
        }

        if (valid) {
            ans++;
        }
    }

    cout << ans << "\n";
    return 0;
}