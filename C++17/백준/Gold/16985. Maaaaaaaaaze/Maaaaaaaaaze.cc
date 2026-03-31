#include <iostream>
#include <vector>
#include <queue>
#include <tuple>
#include <algorithm>

using namespace std;

int original_board[5][5][5];
int board[5][4][5][5]; // [판 번호][회전 상태][y][x]
int dz[] = {1, -1, 0, 0, 0, 0};
int dy[] = {0, 0, 1, -1, 0, 0};
int dx[] = {0, 0, 0, 0, 1, -1};
int min_ans = 999999;

// 각 판의 4가지 회전 상태를 미리 계산해두는 함수
void rotate_board(int idx) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            board[idx][0][i][j] = original_board[idx][i][j];
        }
    }
    for (int r = 1; r < 4; r++) {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                board[idx][r][i][j] = board[idx][r-1][5-1-j][i]; // 90도 시계방향 회전
            }
        }
    }
}

// 3차원 맨해튼 거리 (도착지 4, 4, 4 기준)
int heuristic(int z, int y, int x) {
    return (4 - z) + (4 - y) + (4 - x);
}

void solve() {
    vector<int> order = {0, 1, 2, 3, 4};
    
    // 1. 판을 쌓는 순서 (순열)
    do {
        // 2. 각 판의 회전 상태 (0~1023)
        for (int rot = 0; rot < 1024; rot++) {
            int r[5];
            int temp = rot;
            for (int i = 0; i < 5; i++) {
                r[i] = temp % 4;
                temp /= 4;
            }

            // 입구(0,0,0)와 출구(4,4,4)가 막혀있으면 진행 불가
            if (board[order[0]][r[0]][0][0] == 0 || board[order[4]][r[4]][4][4] == 0) continue;

            // 3. A* 알고리즘 세팅
            // 튜플 구조: {예상 총 비용(f), 지금까지 이동 거리(g), z, y, x}
            priority_queue<tuple<int, int, int, int, int>, 
                           vector<tuple<int, int, int, int, int>>, 
                           greater<tuple<int, int, int, int, int>>> pq;
                           
            int dist[5][5][5];
            for(int i=0; i<5; ++i) 
                for(int j=0; j<5; ++j) 
                    for(int k=0; k<5; ++k) dist[i][j][k] = 999999;

            pq.push({12, 0, 0, 0, 0}); // f=12, g=0 시작
            dist[0][0][0] = 0;

            while (!pq.empty()) {
                int f, g, z, y, x;
                tie(f, g, z, y, x) = pq.top();
                pq.pop();

                // 목적지 도착
                if (z == 4 && y == 4 && x == 4) {
                    min_ans = min(min_ans, g);
                    // 이론상 최단거리면 즉시 종료 (강력한 최적화)
                    if (min_ans == 12) {
                        cout << 12 << "\n";
                        exit(0);
                    }
                    break; 
                }

                // 이미 더 짧은 경로로 방문한 적 있다면 스킵
                if (dist[z][y][x] < g) continue;

                // 6방향 탐색
                for (int dir = 0; dir < 6; dir++) {
                    int nz = z + dz[dir];
                    int ny = y + dy[dir];
                    int nx = x + dx[dir];

                    // 맵 밖으로 나가거나 벽(0)인 경우
                    if (nz < 0 || nz >= 5 || ny < 0 || ny >= 5 || nx < 0 || nx >= 5) continue;
                    if (board[order[nz]][r[nz]][ny][nx] == 0) continue;

                    int h = heuristic(nz, ny, nx);
                    int nf = g + 1 + h;

                    // A* 가지치기: 새로운 예상 비용이 이미 찾은 정답보다 크거나 같으면 가볼 필요 없음
                    if (nf >= min_ans) continue;

                    // 더 짧은 경로를 발견했을 때만 큐에 추가
                    if (g + 1 < dist[nz][ny][nx]) {
                        dist[nz][ny][nx] = g + 1;
                        pq.push({nf, g + 1, nz, ny, nx});
                    }
                }
            }
        }
    } while (next_permutation(order.begin(), order.end()));
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // 입력 받기
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                cin >> original_board[i][j][k];
            }
        }
        // 입력받은 판의 회전 상태 미리 저장
        rotate_board(i);
    }

    solve();

    // 탈출 불가능할 경우 -1 출력
    if (min_ans == 999999) cout << -1 << "\n";
    else cout << min_ans << "\n";

    return 0;
}