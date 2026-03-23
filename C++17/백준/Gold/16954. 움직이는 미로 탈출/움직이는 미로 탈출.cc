#include <iostream>
#include <queue>
#include <string>
#include <vector>
using namespace std;

struct State {
    int r, c, t;
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    // 체스판 입력: 첫 줄은 체스판의 최상단 행, 마지막 줄은 최하단 행
    vector<string> board(8);
    for (int i = 0; i < 8; i++) {
        cin >> board[i];
    }
    
    // 8방향 + 제자리에 머무르는 경우
    int dr[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
    int dc[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    
    // visited[시간][행][열]
    bool visited[9][8][8] = {false, };
    queue<State> q;
    
    // 시작: 가장 왼쪽 아랫칸 (7,0)에서 t=0
    q.push({7,0,0});
    visited[0][7][0] = true;
    
    bool possible = false;
    
    while(!q.empty()){
        State cur = q.front();
        q.pop();
        
        int r = cur.r, c = cur.c, t = cur.t;
        
        // t가 커지면 벽은 모두 사라진다고 볼 수 있음
        if(r == 0 && c == 7) { // 목표지점 도달
            possible = true;
            break;
        }
        
        // 현재 시각의 상태에서 벽이 이미 도달했으면 더 이상 안전하지 않음.
        // (단, t>=8인 경우는 모든 벽이 사라짐)
        if(t < 8 && r - t >= 0 && board[r - t][c] == '#') continue;
        
        // 9가지 이동: 인접 8방향 + 제자리
        for (int i = 0; i < 9; i++){
            int nr = r + dr[i];
            int nc = c + dc[i];
            int nt = t + 1;
            
            // 보드 범위 내 체크
            if(nr < 0 || nr >= 8 || nc < 0 || nc >= 8) continue;
            
            // 현재 시간 t에서 이동할 위치에 벽이 있는지 확인
            if(t < 8 && nr - t >= 0 && board[nr - t][nc] == '#') continue;
            
            // 이동 후, 벽이 내려오는 시간 t+1에서 해당 위치에 벽이 오는지 확인
            if(nt < 8 && nr - nt >= 0 && board[nr - nt][nc] == '#') continue;
            
            if(nt > 8) nt = 8; // 8 이상의 시간은 벽이 모두 사라진 것으로 간주
            
            if(!visited[nt][nr][nc]){
                visited[nt][nr][nc] = true;
                q.push({nr, nc, nt});
            }
        }
    }
    
    cout << (possible ? 1 : 0);
    return 0;
}