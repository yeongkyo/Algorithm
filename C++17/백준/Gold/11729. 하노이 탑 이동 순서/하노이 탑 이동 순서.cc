#include <iostream>
#include <cmath>

using namespace std;

// 하노이 탑 재귀 함수
// n: 원판 개수, from: 출발지, tmp: 보조, to: 목적지
void hanoi(int n, int from, int tmp, int to) {
    if (n == 1) {
        // 원판이 1개면 바로 목적지로 이동
        cout << from << " " << to << "\n";
        return;
    }

    // 1. n-1개를 출발지에서 보조 기둥으로 이동
    hanoi(n - 1, from, to, tmp);

    // 2. 가장 큰 원판 1개를 목적지로 이동
    cout << from << " " << to << "\n";

    // 3. 보조 기둥에 있던 n-1개를 목적지로 이동
    hanoi(n - 1, tmp, from, to);
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    // 총 이동 횟수 K = 2^N - 1
    // (1 << N)은 2의 N제곱을 의미하는 비트 연산입니다.
    cout << (1 << N) - 1 << "\n";

    // 하노이 함수 호출
    hanoi(N, 1, 2, 3);

    return 0;
}