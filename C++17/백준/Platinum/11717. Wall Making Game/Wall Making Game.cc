#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cstring>

using namespace std;

int H, W;
vector<string> board;
int memo[25][25][25][25];

int solve(int r1, int r2, int c1, int c2) {
    if (r1 > r2 || c1 > c2) return 0;
    if (memo[r1][r2][c1][c2] != -1) return memo[r1][r2][c1][c2];

    set<int> mex_set;
    for (int i = r1; i <= r2; ++i) {
        for (int j = c1; j <= c2; ++j) {
            if (board[i][j] == '.') {
                int nim_val = solve(r1, i - 1, c1, j - 1)
                            ^ solve(r1, i - 1, j + 1, c2)
                            ^ solve(i + 1, r2, c1, j - 1)
                            ^ solve(i + 1, r2, j + 1, c2);
                mex_set.insert(nim_val);
            }
        }
    }

    int mex = 0;
    while (mex_set.count(mex)) mex++;

    return memo[r1][r2][c1][c2] = mex;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> H >> W)) return 0;
    board.resize(H);
    for (int i = 0; i < H; ++i) {
        cin >> board[i];
    }

    memset(memo, -1, sizeof(memo));

    if (solve(0, H - 1, 0, W - 1) > 0) {
        cout << "First\n";
    } else {
        cout << "Second\n";
    }

    return 0;
}