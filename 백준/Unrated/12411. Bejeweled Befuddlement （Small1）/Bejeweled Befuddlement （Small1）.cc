#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int simulate(int N, int M, vector<string> board) {
    int total_removed = 0;
    
    while (true) {
        vector<vector<bool>> to_remove(N, vector<bool>(M, false));
        bool found = false;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M - 2; ++j) {
                if (board[i][j] == '.') continue;
                if (board[i][j] == board[i][j+1] && board[i][j] == board[i][j+2]) {
                    found = true;
                    char color = board[i][j];
                    int k = j;
                    while (k < M && board[i][k] == color) {
                        to_remove[i][k] = true;
                        k++;
                    }
                }
            }
        }

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N - 2; ++i) {
                if (board[i][j] == '.') continue;
                if (board[i][j] == board[i+1][j] && board[i][j] == board[i+2][j]) {
                    found = true;
                    char color = board[i][j];
                    int k = i;
                    while (k < N && board[k][j] == color) {
                        to_remove[k][j] = true;
                        k++;
                    }
                }
            }
        }

        if (!found) break;

        int removed_in_step = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                if (to_remove[i][j]) {
                    removed_in_step++;
                    board[i][j] = '.';
                }
            }
        }
        total_removed += removed_in_step;

        for (int j = 0; j < M; ++j) {
            int write_idx = N - 1;
            for (int read_idx = N - 1; read_idx >= 0; --read_idx) {
                if (board[read_idx][j] != '.') {
                    board[write_idx--][j] = board[read_idx][j];
                }
            }
            while (write_idx >= 0) {
                board[write_idx--][j] = '.';
            }
        }
    }
    return total_removed;
}

void solve(int tc) {
    int N, M;
    cin >> N >> M;
    vector<string> board(N);
    for (int i = 0; i < N; ++i) cin >> board[i];

    int max_removed = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M - 1; ++j) {
            swap(board[i][j], board[i][j+1]);
            max_removed = max(max_removed, simulate(N, M, board));
            swap(board[i][j], board[i][j+1]); 
        }
    }

    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < M; ++j) {
            swap(board[i][j], board[i+1][j]);
            max_removed = max(max_removed, simulate(N, M, board));
            swap(board[i][j], board[i+1][j]); 
        }
    }

    cout << "Case #" << tc << ": " << max_removed << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int T;
    cin >> T;
    for (int i = 1; i <= T; ++i) {
        solve(i);
    }
    return 0;
}