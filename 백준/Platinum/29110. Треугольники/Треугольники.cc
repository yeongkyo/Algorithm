#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef long long ll;

struct Point {
    ll x, y;
};

struct Ray {
    ll x, y;
    ll count;
};

int half(Point p) {
    return (p.y > 0 || (p.y == 0 && p.x > 0)) ? 0 : 1;
}

__int128_t cross(Point a, Point b) {
    return (__int128_t)a.x * b.y - (__int128_t)a.y * b.x;
}

__int128_t cross(Ray a, Ray b) {
    return (__int128_t)a.x * b.y - (__int128_t)a.y * b.x;
}

bool cmp(Point a, Point b) {
    if (half(a) != half(b)) return half(a) < half(b);
    __int128_t cp = cross(a, b);
    if (cp != 0) return cp > 0;
    return false; // 외적이 0이면 일직선상에 있음 (나중에 그룹화됨)
}

ll C2(ll n) {
    if (n < 2) return 0;
    return n * (n - 1) / 2;
}

ll C3(ll n) {
    if (n < 3) return 0;
    return n * (n - 1) * (n - 2) / 6;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int n;
    if (!(cin >> n)) return 0;
    
    ll cx, cy;
    cin >> cx >> cy;
    
    vector<Point> pts(n);
    for (int i = 0; i < n; i++) {
        ll x, y;
        cin >> x >> y;
        pts[i] = {x - cx, y - cy}; // 세레브로를 원점으로 평행이동
    }
    
    sort(pts.begin(), pts.end(), cmp);
    
    vector<Ray> rays;
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            rays.push_back({pts[i].x, pts[i].y, 1});
        } else {
            Point a = pts[i-1];
            Point b = pts[i];
            if (half(a) == half(b) && cross(a, b) == 0) {
                rays.back().count++;
            } else {
                rays.push_back({pts[i].x, pts[i].y, 1});
            }
        }
    }
    
    int m = rays.size();
    for (int i = 0; i < m; i++) {
        rays.push_back(rays[i]);
    }
    
    ll bad_triangles = 0;
    int j = 1;
    ll L_sum = 0; // 180도 미만 구간에 속하는 점들의 개수
    
    for (int i = 0; i < m; i++) {
        if (j <= i) {
            j = i + 1;
            L_sum = 0;
        }
        
        while (j < i + m && cross(rays[i], rays[j]) > 0) {
            L_sum += rays[j].count;
            j++;
        }
        
        bad_triangles += C3(rays[i].count + L_sum) - C3(L_sum);
        
        if (j < i + m && cross(rays[i], rays[j]) == 0) {
            if (j < m) {
                ll cA = rays[i].count;
                ll cB = rays[j].count;
                bad_triangles += cA * cB * (n - cA - cB) + C2(cA) * cB + cA * C2(cB);
            }
        }
        
        if (j > i + 1) {
            L_sum -= rays[i + 1].count;
        } else {
            j = i + 2;
            L_sum = 0;
        }
    }
    
    ll total_triangles = C3(n);
    cout << total_triangles - bad_triangles << "\n";
    
    return 0;
}