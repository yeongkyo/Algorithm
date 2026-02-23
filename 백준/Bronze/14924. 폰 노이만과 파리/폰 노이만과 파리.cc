#include <iostream>

using namespace std;

int main() {
    // 입력을 빠르게 받기 위한 설정
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int S, T, D;
    
    // 기차 속도 S, 파리 속도 T, 초기 거리 D 입력
    if (!(cin >> S >> T >> D)) return 0;

    // 1. 두 기차가 만날 때까지 걸리는 시간 계산: 시간 = 거리 / 상대 속도
    // 상대 속도는 두 기차의 속도의 합인 2 * S
    int time = D / (2 * S);

    // 2. 파리가 이동한 거리 계산: 거리 = 속력 * 시간
    int F = T * time;

    // 결과 출력
    cout << F << endl;

    return 0;
}