#include <iostream>

using namespace std;

/**
 * CCW 함수
 * 결과값이 양수면 반시계, 음수면 시계, 0이면 일직선
 */
int ccw(long long x1, long long y1, long long x2, long long y2, long long x3, long long y3) {
    // 외적 공식을 이용한 계산
    // (x1*y2 + x2*y3 + x3*y1) - (y1*x2 + y2*x3 + y3*x1)
    long long temp = (x1 * y2 + x2 * y3 + x3 * y1) - (y1 * x2 + y2 * x3 + y3 * x1);

    if (temp > 0) return 1;      // 반시계 방향
    else if (temp < 0) return -1; // 시계 방향
    else return 0;               // 일직선
}

int main() {
    // 입력 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long x1, y1, x2, y2, x3, y3;

    // 세 점의 좌표 입력
    cin >> x1 >> y1;
    cin >> x2 >> y2;
    cin >> x3 >> y3;

    // 결과 출력
    cout << ccw(x1, y1, x2, y2, x3, y3) << "\n";

    return 0;
}