#include <iostream>
#include <vector>
#include <string>

using namespace std;

int L, W, N, S;
char grid[105][105];

bool solve(int r1, int c1, int r2, int c2, int num, int char_idx) {
    int h = r2 - r1 + 1;
    int w = c2 - c1 + 1;

    if (num == 1) {
        for (int i = r1; i <= r2; ++i) {
            for (int j = c1; j <= c2; ++j) {
                grid[i][j] = 'A' + char_idx;
            }
        }
        return true;
    }

    for (int k = 1; k < num; ++k) {
        long long target_area = (long long)k * S;
        if (target_area % w == 0) {
            int split_h = target_area / w;
            if (split_h > 0 && split_h < h) {
                if (solve(r1, c1, r1 + split_h - 1, c2, k, char_idx) &&
                    solve(r1 + split_h, c1, r2, c2, num - k, char_idx + k)) {
                    return true;
                }
            }
        }
    }

    for (int k = 1; k < num; ++k) {
        long long target_area = (long long)k * S;
        if (target_area % h == 0) {
            int split_w = target_area / h;
            if (split_w > 0 && split_w < w) {
                if (solve(r1, c1, r2, c1 + split_w - 1, k, char_idx) &&
                    solve(r1, c1 + split_w, r2, c2, num - k, char_idx + k)) {
                    return true;
                }
            }
        }
    }

    return false;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> L >> W >> N)) return 0;

    if ((L * W) % N != 0) {
        cout << "impossible" << endl;
        return 0;
    }

    S = (L * W) / N; // 각 직사각형의 고정 면적

    if (solve(0, 0, L - 1, W - 1, N, 0)) {
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < W; ++j) {
                cout << grid[i][j];
            }
            cout << "\n";
        }
    } else {
        cout << "impossible" << endl;
    }

    return 0;
}