#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

struct Point2D {
    double x, y;
    bool operator<(const Point2D& other) const {
        if (abs(x - other.x) > 1e-9) return x < other.x;
        return y < other.y;
    }
};

double cross_product(const Point2D& O, const Point2D& A, const Point2D& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

double convex_hull_area(vector<Point2D>& pts) {
    int n = pts.size(), k = 0;
    if (n < 3) return 0.0;
    vector<Point2D> hull(2 * n);
    sort(pts.begin(), pts.end());

    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross_product(hull[k - 2], hull[k - 1], pts[i]) <= 1e-9) k--;
        hull[k++] = pts[i];
    }

    for (int i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross_product(hull[k - 2], hull[k - 1], pts[i]) <= 1e-9) k--;
        hull[k++] = pts[i];
    }

    hull.resize(k - 1);

    double area = 0.0;
    for (size_t i = 0; i < hull.size(); ++i) {
        size_t j = (i + 1) % hull.size();
        area += hull[i].x * hull[j].y - hull[i].y * hull[j].x;
    }
    return abs(area) / 2.0;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long x1, y1, z1, x2, y2, z2;
    if (!(cin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2)) return 0;

    long long lx, ly, lz;
    cin >> lx >> ly >> lz;

    long long Xmin = min(x1, x2), Xmax = max(x1, x2);
    long long Ymin = min(y1, y2), Ymax = max(y1, y2);
    long long Zmin = min(z1, z2), Zmax = max(z1, z2);

    if (lz <= Zmin) {
        cout << 0 << "\n";
        return 0;
    }

    if (lz <= Zmax && Zmin < lz) {
        if (Xmin < Xmax && Ymin < Ymax) {
            cout << -1 << "\n";
        } else if (Xmin == Xmax && Ymin < Ymax) {
            if (lx == Xmin) cout << 0 << "\n";
            else cout << -1 << "\n";
        } else if (Xmin < Xmax && Ymin == Ymax) {
            if (ly == Ymin) cout << 0 << "\n";
            else cout << -1 << "\n";
        } else {
            cout << 0 << "\n";
        }
        return 0;
    }

    vector<Point2D> pts;
    long long Xs[] = {Xmin, Xmax};
    long long Ys[] = {Ymin, Ymax};
    long long Zs[] = {Zmin, Zmax};

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                double vx = Xs[i];
                double vy = Ys[j];
                double vz = Zs[k];

                double t = (double)lz / (lz - vz);
                double px = lx + t * (vx - lx);
                double py = ly + t * (vy - ly);
                pts.push_back({px, py});
            }
        }
    }

    double area = convex_hull_area(pts);

    if (area < 1e-9) {
        cout << 0 << "\n";
    } else {
        cout << fixed << setprecision(9) << area << "\n";
    }

    return 0;
}