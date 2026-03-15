#include <iostream>
#include <vector>

using namespace std;

const int MOD = 1000000007;

int dp[2][9][10][2][1005];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, A;
    if (!(cin >> N >> A)) return 0;

    // dp[cur][sum_D][last_diff][diff_state][num_blocks]
    // sum_D goes up to 8
    // last_diff goes up to 8, we use 9 to represent the initial state (-1)
    // diff_state is 0 or 1
    // num_blocks is up to A
    
    dp[0][0][9][0][1] = 1;

    int current_max_a = 1;
    for (int i = 0; i < N - 1; ++i) {
        int cur = i % 2;
        int nxt = (i + 1) % 2;

        for (int s = 0; s <= 8; ++s) {
            for (int x = 0; x <= 9; ++x) {
                for (int diff = 0; diff <= 1; ++diff) {
                    for (int a = 1; a <= current_max_a + 1 && a <= A; ++a) {
                        dp[nxt][s][x][diff][a] = 0;
                    }
                }
            }
        }

        for (int s = 0; s <= 8; ++s) {
            for (int x = 0; x <= 9; ++x) {
                for (int diff = 0; diff <= 1; ++diff) {
                    for (int y = 0; y <= 8 - s; ++y) {
                        int nxt_s = s + y;
                        int nxt_diff;
                        int delta_a;
                        
                        if (y == x) {
                            nxt_diff = 1;
                            delta_a = 0;
                        } else {
                            if (diff == 1) {
                                nxt_diff = 0;
                                delta_a = 1;
                            } else {
                                nxt_diff = 1;
                                delta_a = 0;
                            }
                        }

                        int* curr_dp = dp[cur][s][x][diff];
                        int* nxt_dp = dp[nxt][nxt_s][y][nxt_diff];

                        int limit = current_max_a;
                        if (limit > A - delta_a) limit = A - delta_a;

                        for (int a = 1; a <= limit; ++a) {
                            int val = curr_dp[a];
                            if (val) {
                                nxt_dp[a + delta_a] += val;
                                if (nxt_dp[a + delta_a] >= MOD) {
                                    nxt_dp[a + delta_a] -= MOD;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (current_max_a < A) current_max_a++;
    }

    int final_cur = (N - 1) % 2;
    if (N == 1) final_cur = 0;

    long long ans = 0;
    for (int s = 0; s <= 8; ++s) {
        for (int x = 0; x <= 9; ++x) {
            for (int diff = 0; diff <= 1; ++diff) {
                long long ways = dp[final_cur][s][x][diff][A];
                if (ways > 0) {
                    ans = (ans + ways * (9 - s)) % MOD;
                }
            }
        }
    }

    cout << ans << "\n";

    return 0;
}