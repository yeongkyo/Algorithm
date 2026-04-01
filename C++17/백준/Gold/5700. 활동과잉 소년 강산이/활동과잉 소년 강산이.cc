#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// 작업 정보를 담는 구조체 (실수 좌표 대응)
struct Task {
    double s, f;
};

// 시작 시간 기준 오름차순, 같으면 끝 시간 오름차순 정렬
bool compareTasks(const Task& a, const Task& b) {
    if (a.s != b.s) return a.s < b.s;
    return a.f < b.f;
}

// dp[이전_작업_인덱스][현재_작업_인덱스]
long long dp[105][105];
const int MOD = 100000000;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    double M;
    int N;
    // 0 0 입력 시 종료
    while (cin >> M >> N && (M != 0 || N != 0)) {
        vector<Task> tasks(N);
        for (int i = 0; i < N; ++i) {
            cin >> tasks[i].s >> tasks[i].f;
        }

        sort(tasks.begin(), tasks.end(), compareTasks);

        // DP 테이블 초기화 (N은 이전 작업이 없는 경우를 위한 더미 인덱스)
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                dp[i][j] = 0;
            }
        }

        // 시작 시각이 0인 작업들로부터 시작
        for (int i = 0; i < N; ++i) {
            if (tasks[i].s == 0) {
                dp[N][i] = 1;
            }
        }

        // DP 전이
        for (int j = 0; j < N; ++j) { // 현재 마지막 작업
            for (int i = 0; i <= N; ++i) { // 그 직전 작업
                if (dp[i][j] == 0) continue;

                for (int k = j + 1; k < N; ++k) { // 새로 추가할 작업
                    // 1. 시작 시간은 반드시 엄격하게 증가해야 함 (중복 방지)
                    if (tasks[k].s <= tasks[j].s) continue;
                    // 2. 연속성: 현재 작업이 끝난 지점 이전에 다음 작업이 시작되어야 함
                    if (tasks[k].s > tasks[j].f) continue;
                    // 3. 최소성: 직전 직전 작업(i)과 새 작업(k)이 현재 작업(j)을 대체할 수 없어야 함
                    if (i != N && tasks[k].s <= tasks[i].f) continue;

                    dp[j][k] = (dp[j][k] + dp[i][j]) % MOD;
                }
            }
        }

        long long ans = 0;
        // 전체 범위를 덮고 최소성을 유지하는 경우만 합산
        for (int j = 0; j < N; ++j) {
            if (tasks[j].f == M) {
                for (int i = 0; i <= N; ++i) {
                    if (dp[i][j] == 0) continue;
                    // 단일 작업인 경우 시작이 0이고 끝이 M이면 최소 부분집합
                    if (i == N) {
                        if (tasks[j].s == 0) ans = (ans + dp[i][j]) % MOD;
                    } 
                    // 두 개 이상의 작업인 경우 직전 작업이 M에 도달하지 않았어야 현재 작업이 필수적임
                    else {
                        if (tasks[i].f < M) ans = (ans + dp[i][j]) % MOD;
                    }
                }
            }
        }
        cout << ans << "\n";
    }

    return 0;
}