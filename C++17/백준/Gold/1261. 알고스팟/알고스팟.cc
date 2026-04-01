#include <iostream>
#include <vector>
#include <string>
#include <deque>

using namespace std;

// 이동 방향 (상, 하, 좌, 우)
int dx[] = {0, 0, 1, -1};
int dy[] = {1, -1, 0, 0};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int M, N;
    if (!(cin >> M >> N)) return 0;

    vector<string> maze(N);
    for (int i = 0; i < N; i++) {
        cin >> maze[i];
    }

    // dist[r][c]는 (r, c)까지 오기 위해 부순 벽의 최소 개수
    // 초기값은 충분히 큰 값(-1)으로 설정
    vector<vector<int>> dist(N, vector<int>(M, -1));

    deque<pair<int, int>> dq;

    // 시작점 (0, 0) 설정
    dq.push_back({0, 0});
    dist[0][0] = 0;

    while (!dq.empty()) {
        int x = dq.front().first;
        int y = dq.front().second;
        dq.pop_front();

        // 목적지에 도달하면 결과 출력 후 종료
        if (x == N - 1 && y == M - 1) {
            cout << dist[x][y] << "\n";
            return 0;
        }

        for (int i = 0; i < 4; i++) {
            int nx = x + dx[i];
            int ny = y + dy[i];

            // 미로 범위 체크 및 방문 여부 확인
            if (nx >= 0 && nx < N && ny >= 0 && ny < M && dist[nx][ny] == -1) {
                if (maze[nx][ny] == '0') {
                    // 빈 방일 경우: 가중치 0, 덱의 앞에 추가
                    dist[nx][ny] = dist[x][y];
                    dq.push_front({nx, ny});
                } else {
                    // 벽일 경우: 가중치 1, 덱의 뒤에 추가
                    dist[nx][ny] = dist[x][y] + 1;
                    dq.push_back({nx, ny});
                }
            }
        }
    }

    return 0;
}