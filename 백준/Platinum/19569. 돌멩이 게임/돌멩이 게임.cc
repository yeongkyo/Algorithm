#include <iostream>
#include <vector>

using namespace std;

const int MAXN = 100000;
const int MAXK = 450;
bool W[MAXN + 1][MAXK + 1];

int main() {
    // 입출력 속도 향상을 위한 설정 (인터랙티브 문제이므로 endl을 사용해 버퍼를 플러시함)
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    // DP 테이블 채우기
    for (int n = 1; n <= N; ++n) {
        for (int k = 1; k <= MAXK; ++k) {
            // 내가 남은 돌을 다 가져갈 수 있는 경우
            if (k + 1 >= n) {
                W[n][k] = true;
            } else {
                W[n][k] = false;
                // 1개부터 k+1개까지 가져가보면서, 상대방을 지게 만들 수 있는 수가 있는지 확인
                for (int t = 1; t <= k + 1; ++t) {
                    if (!W[n - t][t]) {
                        W[n][k] = true;
                        break;
                    }
                }
            }
        }
    }

    // 첫 번째 턴에는 무조건 1개를 가져가야 하므로, 상대방의 시작 상태는 W[N-1][1] 이 됨.
    // 만약 상대방의 상태가 True 라면 상대가 이기는 게임이므로 나는 패배함.
    if (W[N - 1][1]) {
        cout << "NO" << endl;
        return 0;
    }

    // 내가 이길 수 있는 경우
    cout << "YES" << endl;
    cout << 1 << endl; // 첫 수는 무조건 1

    int rem = N - 1;
    int prev_k = 1;

    // 게임 진행
    while (rem > 0) {
        int x; // muse(상대방)가 가져간 돌의 개수
        if (!(cin >> x)) break;
        
        rem -= x;
        if (rem <= 0) break; // 상대방이 마지막 돌을 가져간 경우 (정상적인 플레이라면 발생하지 않음)
        
        prev_k = x;
        int my_move = -1;

        // 내가 남은 돌을 다 가져갈 수 있다면 다 가져가서 승리
        if (prev_k + 1 >= rem) {
            my_move = rem;
        } else {
            // 상대방을 패배 상태(False)로 만드는 최적의 수 찾기
            for (int t = 1; t <= prev_k + 1; ++t) {
                if (!W[rem - t][t]) {
                    my_move = t;
                    break;
                }
            }
        }

        cout << my_move << endl; // 수 출력 (endl이 버퍼를 자동으로 비워줌)
        rem -= my_move;
        prev_k = my_move;
    }

    return 0;
}