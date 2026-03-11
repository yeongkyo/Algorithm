#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

typedef long long ll;

struct Point {
    ll x, y;
};

int main() {
    // 입출력 속도 향상
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    ll l;
    if (!(cin >> n >> l)) return 0;

    vector<Point> p(n);
    for (int i = 0; i < n; ++i) {
        cin >> p[i].x >> p[i].y;
    }

    ll max_2a = 0;
    ll limit_2a = 2 * l; // 소수점 오차 방지를 위해 면적에 2를 곱해서 정수로 계산

    // 메모리 재할당 방지를 위해 밖에서 선언
    vector<ll> side1, side2;
    side1.reserve(n);
    side2.reserve(n);

    // 1. 모든 점 쌍 (i, j)를 사각형의 대각선으로 고정
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            side1.clear();
            side2.clear();

            // 2. 대각선을 기준으로 나머지 점들을 두 그룹으로 나눔 (외적 사용)
            for (int k = 0; k < n; ++k) {
                if (k == i || k == j) continue;
                
                ll cp = (p[j].x - p[i].x) * (p[k].y - p[i].y) - (p[j].y - p[i].y) * (p[k].x - p[i].x);
                
                if (cp > 0) side1.push_back(cp);
                else if (cp < 0) side2.push_back(-cp); // 양수로 저장
            }

            if (side1.empty() || side2.empty()) continue;

            sort(side1.begin(), side1.end());
            sort(side2.begin(), side2.end());

            int ptr2 = (int)side2.size() - 1;
            for (int ptr1 = 0; ptr1 < (int)side1.size(); ++ptr1) {
                while (ptr2 >= 0 && side1[ptr1] + side2[ptr2] > limit_2a) {
                    ptr2--;
                }
                if (ptr2 >= 0) {
                    ll current_2a = side1[ptr1] + side2[ptr2];
                    if (current_2a > max_2a) {
                        max_2a = current_2a;
                    }
                }
            }
        }
    }

    if (max_2a == 0) {
        cout << "0.00\n";
    } else {
        cout << fixed << setprecision(2) << (double)max_2a / 2.0 << "\n";
    }

    return 0;
}