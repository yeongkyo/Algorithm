#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// 2차원 좌표를 나타내는 구조체
struct Point {
    double x, y;
};

// 특정 각도(theta)로 회전된 격자계에서의 맨해튼 거리 계산
double get_dist(Point p1, Point p2, double theta) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    
    // 회전 변환 행렬 적용
    double X = dx * cos(theta) + dy * sin(theta);
    double Y = -dx * sin(theta) + dy * cos(theta);
    
    return abs(X) + abs(Y);
}

int main() {
    // 입출력 속도 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<Point> pts(n);
    for (int i = 0; i < n; ++i) {
        cin >> pts[i].x >> pts[i].y;
    }

    // 최적의 격자 방향이 될 수 있는 후보 각도 수집
    vector<double> angles;
    angles.push_back(0.0); // 기본 각도
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dy = pts[j].y - pts[i].y;
            double dx = pts[j].x - pts[i].x;
            angles.push_back(atan2(dy, dx)); // 두 점을 잇는 선분의 각도
        }
    }

    double min_total_dist = 1e18; // 무한대 값으로 초기화

    // 각 후보 각도별로 TSP(외판원 순회 문제) 실행
    for (double theta : angles) {
        
        // 미리 현재 각도에서의 모든 점 쌍의 거리를 계산해 둠 (최적화)
        vector<vector<double>> dist(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                dist[i][j] = get_dist(pts[i], pts[j], theta);
            }
        }

        // DP 배열 초기화: dp[방문한 정점 비트마스크][현재 위치] = 최소 거리
        vector<vector<double>> dp(1 << n, vector<double>(n, 1e18));
        
        // 시작점 초기화 (어느 점에서든 시작할 수 있음)
        for (int i = 0; i < n; ++i) {
            dp[1 << i][i] = 0.0;
        }

        // 비트마스크 DP 진행
        for (int mask = 1; mask < (1 << n); ++mask) {
            for (int i = 0; i < n; ++i) {
                // 현재 정점 i를 방문하지 않은 상태라면 스킵
                if (!(mask & (1 << i))) continue;
                if (dp[mask][i] == 1e18) continue;

                for (int j = 0; j < n; ++j) {
                    // 다음 정점 j를 이미 방문했다면 스킵
                    if (mask & (1 << j)) continue;
                    
                    double next_cost = dp[mask][i] + dist[i][j];
                    if (next_cost < dp[mask | (1 << j)][j]) {
                        dp[mask | (1 << j)][j] = next_cost;
                    }
                }
            }
        }

        // 모든 점을 방문한 상태의 최솟값 갱신
        for (int i = 0; i < n; ++i) {
            min_total_dist = min(min_total_dist, dp[(1 << n) - 1][i]);
        }
    }

    // 결과 출력 (오차 허용 범위 내인 소수점 11자리까지)
    cout << fixed << setprecision(11) << min_total_dist << "\n";

    return 0;
}