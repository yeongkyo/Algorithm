#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int N;
    if (!(cin >> N)) return 0;

    vector<vector<int>> board(N, vector<int>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cin >> board[i][j];
        }
    }

    int max_sum = -127 * 100 * 100;

    for (int i = 0; i < N; ++i) {
        vector<int> temp(N, 0);
        for (int j = i; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                temp[k] += board[j][k];
            }

            int current_max = temp[0];
            int current_running = temp[0];
            for (int k = 1; k < N; ++k) {
                current_running = max(temp[k], current_running + temp[k]);
                current_max = max(current_max, current_running);
            }
            max_sum = max(max_sum, current_max);
        }
    }

    cout << max_sum << endl;

    return 0;
}