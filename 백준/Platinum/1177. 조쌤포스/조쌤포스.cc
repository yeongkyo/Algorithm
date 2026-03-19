#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct Event {
    double t;
    int type;
    bool operator<(const Event& other) const {
        if (abs(t - other.t) > 1e-11) {
            return t < other.t;
        }
        return type > other.type;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    double r, bx, by, bvx, bvy;
    if (!(cin >> n >> r >> bx >> by >> bvx >> bvy)) return 0;

    vector<Event> events;
    int base_count = 0;
    double Reff = r + 0.0001;
    double Reff_sq = Reff * Reff;

    for (int i = 0; i < n; i++) {
        double x, y, vx, vy;
        cin >> x >> y >> vx >> vy;

        double dx = x - bx;
        double dy = y - by;
        double dvx = vx - bvx;
        double dvy = vy - bvy;

        double a = dvx * dvx + dvy * dvy;
        double b = 2.0 * (dx * dvx + dy * dvy);
        double c = dx * dx + dy * dy - Reff_sq;

        if (a == 0.0) {
            if (c <= 0.0) {
                base_count++;
            }
        } else {
            double D = b * b - 4.0 * a * c;
            if (D >= 0.0) {
                double sqrt_D = sqrt(D);
                double t1 = (-b - sqrt_D) / (2.0 * a);
                double t2 = (-b + sqrt_D) / (2.0 * a);

                if (t2 >= 0.0) {
                    double t_start = max(0.0, t1);
                    double t_end = t2;
                    events.push_back({t_start, 1});
                    events.push_back({t_end, -1});
                }
            }
        }
    }

    sort(events.begin(), events.end());

    int current = base_count;
    int max_caught = base_count;

    for (const auto& ev : events) {
        current += ev.type;
        if (current > max_caught) {
            max_caught = current;
        }
    }

    cout << max_caught << "\n";

    return 0;
}