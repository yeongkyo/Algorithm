#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<string> board(N);
    for (int i = 0; i < N; ++i) {
        cin >> board[i];
    }

    int min_t = N * N;

    for (int step = 0; step < (1 << N); ++step) {
        int sum = 0;
        for (int j = 0; j < N; ++j) {
            int t_cnt = 0;
            for (int i = 0; i < N; ++i) {
                char current = board[i][j];
                if ((step & (1 << i)) != 0) {
                    current = (current == 'H') ? 'T' : 'H';
                }
                if (current == 'T') {
                    t_cnt++;
                }
            }
            sum += min(t_cnt, N - t_cnt);
        }
        min_t = min(min_t, sum);
    }

    cout << min_t << "\n";

    return 0;
}