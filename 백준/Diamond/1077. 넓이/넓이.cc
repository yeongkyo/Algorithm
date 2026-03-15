#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

struct Point {
    double x, y;
};

// 벡터의 외적을 계산하는 함수
double cross_product(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// 점이 선분 위에 있는지 확인하는 함수
bool is_point_on_segment(Point p, Point a, Point b) {
    return cross_product(a, b, p) == 0 && (p.x >= min(a.x, b.x) && p.x <= max(a.x, b.x)) && (p.y >= min(a.y, b.y) && p.y <= max(a.y, b.y));
}

// 두 선분의 교점을 찾는 함수
Point intersection(Point a, Point b, Point c, Point d) {
    double det = (b.x - a.x) * (d.y - c.y) - (b.y - a.y) * (d.x - c.x);
    if (det == 0) {
        return {NAN, NAN}; // 평행한 경우
    }
    double t = ((c.x - a.x) * (d.y - c.y) - (c.y - a.y) * (d.x - c.x)) / det;
    double u = -((a.x - c.x) * (b.y - a.y) - (a.y - c.y) * (b.x - a.x)) / det;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
    }
    return {NAN, NAN}; // 교점이 없는 경우
}

// 다각형의 넓이를 계산하는 함수
double polygon_area(const vector<Point>& polygon) {
    double area = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i) {
        area += (polygon[i].x * polygon[(i + 1) % polygon.size()].y - polygon[(i + 1) % polygon.size()].x * polygon[i].y);
    }
    return abs(area) / 2.0;
}

int main() {
    int n, m;
    cin >> n >> m;

    vector<Point> polygon1(n);
    for (int i = 0; i < n; ++i) {
        cin >> polygon1[i].x >> polygon1[i].y;
    }

    vector<Point> polygon2(m);
    for (int i = 0; i < m; ++i) {
        cin >> polygon2[i].x >> polygon2[i].y;
    }

    vector<Point> intersection_polygon;

    // 다각형1의 각 점이 다각형2 내부에 있는지 확인
    for (int i = 0; i < n; ++i) {
        bool inside = true;
        for (int j = 0; j < m; ++j) {
            if (cross_product(polygon2[j], polygon2[(j + 1) % m], polygon1[i]) < 0) {
                inside = false;
                break;
            }
        }
        if (inside) {
            intersection_polygon.push_back(polygon1[i]);
        }
    }

    // 다각형2의 각 점이 다각형1 내부에 있는지 확인
    for (int i = 0; i < m; ++i) {
        bool inside = true;
        for (int j = 0; j < n; ++j) {
            if (cross_product(polygon1[j], polygon1[(j + 1) % n], polygon2[i]) < 0) {
                inside = false;
                break;
            }
        }
        if (inside) {
            intersection_polygon.push_back(polygon2[i]);
        }
    }

    // 다각형1의 각 변과 다각형2의 각 변의 교점을 찾음
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            Point p = intersection(polygon1[i], polygon1[(i + 1) % n], polygon2[j], polygon2[(j + 1) % m]);
            if (!isnan(p.x)) {
                intersection_polygon.push_back(p);
            }
        }
    }

    // 교차 다각형의 점들을 정렬하여 볼록 다각형을 만듦
    if (!intersection_polygon.empty()) {
        Point center = {0, 0};
        for (const auto& p : intersection_polygon) {
            center.x += p.x;
            center.y += p.y;
        }
        center.x /= intersection_polygon.size();
        center.y /= intersection_polygon.size();

        sort(intersection_polygon.begin(), intersection_polygon.end(), [&](Point a, Point b) {
            return atan2(a.y - center.y, a.x - center.x) < atan2(b.y - center.y, b.x - center.x);
        });
        
        // 중복 점 제거
        auto last = unique(intersection_polygon.begin(), intersection_polygon.end(), [](Point a, Point b){
            return abs(a.x - b.x) < 1e-10 && abs(a.y-b.y) < 1e-10;
        });
        intersection_polygon.erase(last, intersection_polygon.end());
    }

    cout << fixed << setprecision(15) << polygon_area(intersection_polygon) << endl;

    return 0;
}