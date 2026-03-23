#include <iostream>
#include <vector>
#include <string>

using namespace std;

int R, C;
string grid[25];

bool pointsDown(char b) { return b == '|' || b == '+' || b == '1' || b == '4'; }
bool pointsUp(char b) { return b == '|' || b == '+' || b == '2' || b == '3'; }
bool pointsRight(char b) { return b == '-' || b == '+' || b == '1' || b == '2'; }
bool pointsLeft(char b) { return b == '-' || b == '+' || b == '3' || b == '4'; }

bool needsConnection(int r, int c) {
    if (grid[r][c] == 'M' || grid[r][c] == 'Z') {
        int dr[] = {-1, 1, 0, 0};
        int dc[] = {0, 0, -1, 1};
        for (int i = 0; i < 4; ++i) {
            int nr = r + dr[i];
            int nc = c + dc[i];
            if (nr < 0 || nr >= R || nc < 0 || nc >= C) continue;
            char b = grid[nr][nc];
            if (b == '.') continue;
            if (nr < r && pointsDown(b)) return false;
            if (nr > r && pointsUp(b)) return false;
            if (nc < c && pointsRight(b)) return false;
            if (nc > c && pointsLeft(b)) return false;
        }
        return true;
    }
    return false;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> R >> C;
    for (int i = 0; i < R; ++i) cin >> grid[i];

    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (grid[r][c] == '.') {
                bool u_p = (r > 0 && pointsDown(grid[r - 1][c]));
                bool d_p = (r < R - 1 && pointsUp(grid[r + 1][c]));
                bool l_p = (c > 0 && pointsRight(grid[r][c - 1]));
                bool r_p = (c < C - 1 && pointsLeft(grid[r][c + 1]));

                bool u_m = (r > 0 && needsConnection(r - 1, c));
                bool d_m = (r < R - 1 && needsConnection(r + 1, c));
                bool l_m = (c > 0 && needsConnection(r, c - 1));
                bool r_m = (c < C - 1 && needsConnection(r, c + 1));

                if (u_p || d_p || l_p || r_p || u_m || d_m || l_m || r_m) {
                    bool u = u_p || u_m;
                    bool d = d_p || d_m;
                    bool l = l_p || l_m;
                    bool rr = r_p || r_m;

                    char res;
                    if (u && d && l && rr) res = '+';
                    else if (u && d) res = '|';
                    else if (l && rr) res = '-';
                    else if (d && rr) res = '1';
                    else if (u && rr) res = '2';
                    else if (u && l) res = '3';
                    else if (d && l) res = '4';

                    cout << r + 1 << " " << c + 1 << " " << res << endl;
                    return 0;
                }
            }
        }
    }
    return 0;
}