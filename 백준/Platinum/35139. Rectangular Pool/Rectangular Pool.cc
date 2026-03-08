#include <iostream>
#include <vector>
#include <string>

using namespace std;

// 전역 변수로 그리드와 크기, 나쁜 패턴의 개수 선언
int N, M;
vector<string> grid;
int bad_count = 0;

// (r, c)를 좌측 상단으로 하는 2x2 영역이 '나쁜 패턴'인지 확인하는 함수
// 나쁜 패턴: 2x2 영역 내에 웅덩이('.')가 정확히 3개인 경우
bool is_bad(int r, int c) {
    // 범위를 벗어나면 false (2x2를 형성할 수 없음)
    if (r < 0 || r >= N - 1 || c < 0 || c >= M - 1) return false;
    
    int puddles = 0;
    if (grid[r][c] == '.') puddles++;
    if (grid[r+1][c] == '.') puddles++;
    if (grid[r][c+1] == '.') puddles++;
    if (grid[r+1][c+1] == '.') puddles++;
    
    return puddles == 3;
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M)) return 0;

    grid.resize(N);
    for (int i = 0; i < N; ++i) {
        cin >> grid[i];
    }

    // 초기 상태에서 나쁜 패턴의 개수를 셈
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < M - 1; ++j) {
            if (is_bad(i, j)) {
                bad_count++;
            }
        }
    }

    int Q;
    cin >> Q;
    while (Q--) {
        int r, c;
        cin >> r >> c;
        // 1-based index를 0-based index로 변환
        --r; --c;

        // (r, c) 셀의 변화는 최대 4개의 2x2 영역에 영향을 줌
        // 해당 영역들의 좌측 상단 좌표: (r-1, c-1), (r-1, c), (r, c-1), (r, c)
        int dr[] = {-1, -1, 0, 0};
        int dc[] = {-1, 0, -1, 0};

        // 1. 변경 전: 영향을 받는 2x2 영역 중 나쁜 패턴이 있다면 카운트 감소
        for (int k = 0; k < 4; ++k) {
            if (is_bad(r + dr[k], c + dc[k])) {
                bad_count--;
            }
        }

        // 2. 셀 상태 변경
        if (grid[r][c] == '.') {
            grid[r][c] = '#';
        } else {
            grid[r][c] = '.';
        }

        // 3. 변경 후: 영향을 받는 2x2 영역 중 나쁜 패턴이 생겼다면 카운트 증가
        for (int k = 0; k < 4; ++k) {
            if (is_bad(r + dr[k], c + dc[k])) {
                bad_count++;
            }
        }

        // 나쁜 패턴이 하나도 없다면 모든 웅덩이는 직사각형임
        if (bad_count == 0) {
            cout << "RECTANGLES\n";
        } else {
            cout << "NO\n";
        }
    }

    return 0;
}