#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int n, m, k;
char board[2000][2000];
int psum[2001][2001];

int solve(char color) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int diff = 0;
            if ((i + j) % 2 == 0) {
                if (board[i][j] != color) diff = 1;
            } else {
                if (board[i][j] == color) diff = 1;
            }
            psum[i + 1][j + 1] = psum[i][j + 1] + psum[i + 1][j] - psum[i][j] + diff;
        }
    }

    int result = 2e9;
    for (int i = k; i <= n; i++) {
        for (int j = k; j <= m; j++) {
            int cnt = psum[i][j] - psum[i - k][j] - psum[i][j - k] + psum[i - k][j - k];
            if (cnt < result) result = cnt;
        }
    }
    return result;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> n >> m >> k;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> board[i][j];
        }
    }

    int res1 = solve('B');
    int res2 = solve('W');

    cout << (res1 < res2 ? res1 : res2) << "\n";

    return 0;
}