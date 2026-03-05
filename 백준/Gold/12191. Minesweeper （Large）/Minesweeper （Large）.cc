#include <iostream>
#include <vector>
#include <string>
#include <queue>

using namespace std;

int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};
int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};

void solve(int tc) {
    int N;
    cin >> N;
    vector<string> grid(N);
    for (int i = 0; i < N; ++i) {
        cin >> grid[i];
    }

    vector<vector<int>> mines(N, vector<int>(N, 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (grid[i][j] == '*') continue; 
            
            int count = 0;
            for (int d = 0; d < 8; ++d) {
                int ni = i + dx[d];
                int nj = j + dy[d];
                if (ni >= 0 && ni < N && nj >= 0 && nj < N && grid[ni][nj] == '*') {
                    count++;
                }
            }
            mines[i][j] = count;
        }
    }

    vector<vector<bool>> visited(N, vector<bool>(N, false));
    int clicks = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (grid[i][j] == '.' && mines[i][j] == 0 && !visited[i][j]) {
                clicks++; 
                queue<pair<int, int>> q;
                q.push({i, j});
                visited[i][j] = true;

                while (!q.empty()) {
                    auto [r, c] = q.front();
                    q.pop();

                    for (int d = 0; d < 8; ++d) {
                        int nr = r + dx[d];
                        int nc = c + dy[d];
                        
                        if (nr >= 0 && nr < N && nc >= 0 && nc < N && !visited[nr][nc] && grid[nr][nc] == '.') {
                            visited[nr][nc] = true; 
                            
                            if (mines[nr][nc] == 0) {
                                q.push({nr, nc});
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (grid[i][j] == '.' && !visited[i][j]) {
                clicks++;
            }
        }
    }

    cout << "Case #" << tc << ": " << clicks << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        for (int i = 1; i <= T; ++i) {
            solve(i);
        }
    }
    return 0;
}