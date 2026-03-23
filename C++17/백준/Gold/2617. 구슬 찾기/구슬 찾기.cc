#include <iostream>
#include <vector>

using namespace std;

bool dist[101][101];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M;
    cin >> N >> M;

    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        dist[u][v] = true; // u가 v보다 무겁다
    }

    // 플로이드-워셜 알고리즘으로 전이 관계 형성
    for (int k = 1; k <= N; ++k) {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                if (dist[i][k] && dist[k][j]) {
                    dist[i][j] = true;
                }
            }
        }
    }

    int mid = (N + 1) / 2;
    int answer = 0;

    for (int i = 1; i <= N; ++i) {
        int heavier = 0;
        int lighter = 0;

        for (int j = 1; j <= N; ++j) {
            if (i == j) continue;
            if (dist[j][i]) heavier++;
            if (dist[i][j]) lighter++; // i가 j보다 무거움
        }

        if (heavier >= mid || lighter >= mid) {
            answer++;
        }
    }

    cout << answer << "\n";

    return 0;
}