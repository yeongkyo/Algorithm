#include <iostream>
#include <vector>

using namespace std;

int N, M, K;
int notebook[42][42]; // 노트북의 상태 (1: 스티커 있음, 0: 빈 칸)

// 스티커를 시계 방향으로 90도 회전시키는 함수
void rotate(vector<vector<int>>& sticker) {
    int R = sticker.size();
    int C = sticker[0].size();
    
    // 회전된 스티커의 크기는 C x R
    vector<vector<int>> new_sticker(C, vector<int>(R));
    
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            // (i, j) -> (j, R-1-i) 로 이동
            new_sticker[j][R - 1 - i] = sticker[i][j];
        }
    }
    sticker = new_sticker; // 원본 스티커 교체
}

// (r, c) 위치에 스티커를 붙일 수 있는지 확인하는 함수
bool can_attach(int r, int c, const vector<vector<int>>& sticker) {
    int R = sticker.size();
    int C = sticker[0].size();
    
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            // 스티커가 있는 부분(1)이 노트북의 이미 붙은 부분(1)과 겹치면 불가
            if (sticker[i][j] == 1 && notebook[r + i][c + j] == 1) {
                return false;
            }
        }
    }
    return true;
}

// (r, c) 위치에 스티커를 실제로 붙이는 함수
void attach(int r, int c, const vector<vector<int>>& sticker) {
    int R = sticker.size();
    int C = sticker[0].size();
    
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) {
            if (sticker[i][j] == 1) {
                notebook[r + i][c + j] = 1;
            }
        }
    }
}

int main() {
    // 입출력 속도 향상
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> N >> M >> K;

    while (K--) {
        int R, C;
        cin >> R >> C;
        
        // 스티커 입력 받기
        vector<vector<int>> sticker(R, vector<int>(C));
        for (int i = 0; i < R; ++i) {
            for (int j = 0; j < C; ++j) {
                cin >> sticker[i][j];
            }
        }

        // 4방향 회전 (0, 90, 180, 270도)
        bool attached = false;
        for (int rot = 0; rot < 4; ++rot) {
            R = sticker.size();
            C = sticker[0].size();

            // 노트북의 모든 위치 탐색 (위쪽 -> 왼쪽 순서)
            // 주의: 스티커 크기(R, C)가 노트북(N, M)보다 클 수 있으므로 범위 체크 중요
            // 스티커가 붙을 수 있는 시작점은 (0,0) ~ (N-R, M-C) 까지임
            for (int r = 0; r <= N - R; ++r) {
                for (int c = 0; c <= M - C; ++c) {
                    if (can_attach(r, c, sticker)) {
                        attach(r, c, sticker);
                        attached = true;
                        break; // 해당 스티커 붙였으면 탐색 종료
                    }
                }
                if (attached) break;
            }
            
            if (attached) break; // 붙였으면 다음 스티커로 이동
            rotate(sticker); // 못 붙였으면 90도 회전
        }
    }

    // 최종적으로 노트북에 붙은 스티커 칸의 수 세기
    int cnt = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            cnt += notebook[i][j];
        }
    }
    
    cout << cnt << "\n";

    return 0;
}