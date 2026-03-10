#include <iostream>
#include <vector>

using namespace std;

bool dp[301][301];
int grid[301][301];

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cin >> grid[i][j];
        }
    }

    dp[0][0] = true;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (grid[i][j] == 0) continue;

            if (i > 0 && dp[i - 1][j]) {
                dp[i][j] = true;
            }
            if (j > 0 && dp[i][j - 1]) {
                dp[i][j] = true;
            }
        }
    }

    if (dp[m - 1][n - 1]) {
        cout << "Yes" << "\n";
    } else {
        cout << "No" << "\n";
    }

    return 0;
}