#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>

using namespace std;

struct Point {
    double x, y;
    Point operator-(const Point& o) const { return {x - o.x, y - o.y}; }
    Point operator+(const Point& o) const { return {x + o.x, y + o.y}; }
};

double dot(Point a, Point b) { return a.x * b.x + a.y * b.y; }
double cross(Point a, Point b) { return a.x * b.y - a.y * b.x; }
double dist(Point a, Point b) { return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)); }

struct Polygon {
    vector<Point> pts;
    double min_x, max_x, min_y, max_y;
};

double distToSegment(Point P, Point A, Point B) {
    Point AB = B - A;
    Point AP = P - A;
    double l2 = dot(AB, AB);
    if (l2 < 1e-9) return dist(P, A);
    double t = dot(AP, AB) / l2;
    t = max(0.0, min(1.0, t));
    Point proj = {A.x + t * AB.x, A.y + t * AB.y};
    return dist(P, proj);
}

bool isStrictlyInside(Point P, const Polygon& poly) {
    if (P.x < poly.min_x - 1e-7 || P.x > poly.max_x + 1e-7 ||
        P.y < poly.min_y - 1e-7 || P.y > poly.max_y + 1e-7) return false;

    for (size_t i = 0; i < poly.pts.size(); i++) {
        Point A = poly.pts[i];
        Point B = poly.pts[(i + 1) % poly.pts.size()];
        if (distToSegment(P, A, B) < 1e-7) return false;
    }

    int count = 0;
    for (size_t i = 0; i < poly.pts.size(); i++) {
        Point A = poly.pts[i];
        Point B = poly.pts[(i + 1) % poly.pts.size()];
        if (A.y > B.y) swap(A, B);
        if (P.y < A.y || P.y >= B.y) continue;
        if (abs(A.y - B.y) < 1e-9) continue;
        double x_int = A.x + (P.y - A.y) * (B.x - A.x) / (B.y - A.y);
        if (x_int > P.x) count++;
    }
    return count % 2 == 1;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    double N, M;
    if (!(cin >> N >> M)) return 0;

    Point S;
    cin >> S.x >> S.y;

    int G;
    cin >> G;

    vector<Polygon> polys(G);
    for (int i = 0; i < G; i++) {
        int K;
        cin >> K;
        double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
        for (int j = 0; j < K; j++) {
            Point p;
            cin >> p.x >> p.y;
            polys[i].pts.push_back(p);
            min_x = min(min_x, p.x);
            max_x = max(max_x, p.x);
            min_y = min(min_y, p.y);
            max_y = max(max_y, p.y);
        }
        polys[i].min_x = min_x;
        polys[i].max_x = max_x;
        polys[i].min_y = min_y;
        polys[i].max_y = max_y;
    }

    vector<Point> candidates;
    candidates.push_back({0, 0});
    candidates.push_back({N, 0});
    candidates.push_back({N, M});
    candidates.push_back({0, M});

    auto addCandidate = [&](Point V) {
        Point D = V - S;
        if (abs(D.x) < 1e-9 && abs(D.y) < 1e-9) return;
        double min_t = 1e12;
        if (D.x > 1e-9) min_t = min(min_t, (N - S.x) / D.x);
        if (D.x < -1e-9) min_t = min(min_t, (0 - S.x) / D.x);
        if (D.y > 1e-9) min_t = min(min_t, (M - S.y) / D.y);
        if (D.y < -1e-9) min_t = min(min_t, (0 - S.y) / D.y);

        if (min_t > 0 && min_t < 1e11) {
            Point E = {S.x + min_t * D.x, S.y + min_t * D.y};
            E.x = max(0.0, min(N, E.x));
            E.y = max(0.0, min(M, E.y));
            candidates.push_back(E);
        }
    };

    for (const auto& poly : polys) {
        for (const auto& V : poly.pts) {
            addCandidate(V);
        }
    }

    double max_dist = -1.0;
    Point best_E = {-1, -1};

    for (Point E : candidates) {
        Point vec_v = E - S;
        double len_v = sqrt(dot(vec_v, vec_v));
        if (len_v < 1e-9) continue;

        double min_sx = min(S.x, E.x), max_sx = max(S.x, E.x);
        double min_sy = min(S.y, E.y), max_sy = max(S.y, E.y);

        vector<double> ts;
        ts.push_back(0.0);
        ts.push_back(1.0);

        for (const auto& poly : polys) {
            if (max_sx < poly.min_x - 1e-7 || min_sx > poly.max_x + 1e-7 ||
                max_sy < poly.min_y - 1e-7 || min_sy > poly.max_y + 1e-7) continue;

            for (const auto& V : poly.pts) {
                Point vec_SV = V - S;
                double cr = cross(vec_v, vec_SV);
                double d_to_line = abs(cr) / len_v;
                if (d_to_line < 1e-7) {
                    double t = dot(vec_SV, vec_v) / (len_v * len_v);
                    if (t >= -1e-9 && t <= 1.0 + 1e-9) {
                        ts.push_back(max(0.0, min(1.0, t)));
                    }
                }
            }

            for (size_t i = 0; i < poly.pts.size(); i++) {
                Point A = poly.pts[i];
                Point B = poly.pts[(i + 1) % poly.pts.size()];
                Point vec_w = B - A;
                double denom = cross(vec_v, vec_w);
                if (abs(denom) > 1e-9) {
                    double t = cross(A - S, vec_w) / denom;
                    double u = cross(A - S, vec_v) / denom;
                    if (t >= -1e-9 && t <= 1.0 + 1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9) {
                        ts.push_back(max(0.0, min(1.0, t)));
                    }
                }
            }
        }

        sort(ts.begin(), ts.end());
        vector<double> unique_ts;
        for (double t : ts) {
            if (unique_ts.empty() || t - unique_ts.back() > 1e-7) {
                unique_ts.push_back(t);
            }
        }

        bool isValid = true;
        for (size_t i = 0; i < unique_ts.size() - 1; i++) {
            double mid_t = (unique_ts[i] + unique_ts[i+1]) / 2.0;
            Point M_pt = {S.x + mid_t * vec_v.x, S.y + mid_t * vec_v.y};
            for (const auto& poly : polys) {
                if (isStrictlyInside(M_pt, poly)) {
                    isValid = false;
                    break;
                }
            }
            if (!isValid) break;
        }

        if (isValid) {
            double d = dist(S, E);
            if (d > max_dist + 1e-7) {
                max_dist = d;
                best_E = E;
            } else if (abs(d - max_dist) <= 1e-7) {
                if (E.x < best_E.x - 1e-7) {
                    best_E = E;
                } else if (abs(E.x - best_E.x) <= 1e-7) {
                    if (E.y < best_E.y - 1e-7) {
                        best_E = E;
                    }
                }
            }
        }
    }

    if (max_dist < 0) {
        cout << "GG\n";
    } else {
        printf("%.3f %.3f\n", abs(best_E.x) < 1e-9 ? 0.0 : best_E.x, abs(best_E.y) < 1e-9 ? 0.0 : best_E.y);
    }

    return 0;
}