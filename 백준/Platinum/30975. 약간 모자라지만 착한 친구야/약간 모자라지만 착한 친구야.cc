#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int INF = 1e9; // 도달할 수 없음을 나타내는 충분히 큰 값

int main() {
    // 입출력 속도 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<int> p(n);
    for (int i = 0; i < n; ++i) {
        cin >> p[i];
        p[i]--; // 0-based index로 변환
    }

    // 0 ~ n-1: 일반 동네, n: 경상국립대학교
    vector<vector<int>> dist(n + 1, vector<int>(n + 1, INF));
    for (int i = 0; i < m; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        u--; v--; // 0-based index로 변환
        dist[u][v] = min(dist[u][v], w); // 동일한 경로의 여러 간선 중 최솟값만 저장
    }

    // dp[방문한 동네 비트마스크][현재 위치한 동네] = 최소 시간
    vector<vector<int>> dp(1 << n, vector<int>(n, INF));

    // 1. 대학교에서 처음 출발하는 경우 초기화
    for (int i = 0; i < n; ++i) {
        // 선행 조건이 없고(자기 자신이 조건), 대학교에서 가는 길이 있을 때만 시작 가능
        if (p[i] == i && dist[n][i] != INF) {
            dp[1 << i][i] = dist[n][i];
        }
    }

    // 2. 비트마스크 DP 탐색
    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int u = 0; u < n; ++u) {
            // 현재 동네 u를 방문하지 않은 상태라면 건너뜀
            if (!(mask & (1 << u))) continue;
            // 도달할 수 없는 상태면 건너뜀
            if (dp[mask][u] == INF) continue;

            for (int v = 0; v < n; ++v) {
                // 다음 동네 v를 이미 방문했다면 건너뜀
                if (mask & (1 << v)) continue;
                // u에서 v로 가는 간선이 없다면 건너뜀
                if (dist[u][v] == INF) continue;

                // 다음 동네 v를 방문하기 위한 선행 조건 검사
                // 선행 동네가 자기 자신이 아니고, 아직 방문하지 않았다면 갈 수 없음
                if (p[v] != v && !(mask & (1 << p[v]))) continue;

                // 조건 만족 시 최솟값 갱신
                dp[mask | (1 << v)][v] = min(dp[mask | (1 << v)][v], dp[mask][u] + dist[u][v]);
            }
        }
    }

    // 3. 모든 동네를 방문한 후 다시 대학교로 돌아가는 최솟값 계산
    int ans = INF;
    int full_mask = (1 << n) - 1;
    
    for (int i = 0; i < n; ++i) {
        // 모든 동네를 방문했고, 마지막 동네 i에서 대학교 n으로 가는 길이 존재할 때
        if (dp[full_mask][i] != INF && dist[i][n] != INF) {
            ans = min(ans, dp[full_mask][i] + dist[i][n]);
        }
    }

    // 결과 출력
    if (ans == INF) cout << -1 << "\n";
    else cout << ans << "\n";

    return 0;
}