#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

int N, M;
long long C[1005][1005];

void print_grid(const vector<long long>& candR, const vector<long long>& candC) {
    vector<vector<long long>> A(N, vector<long long>(M, 0));
    long long sum_C_except_last = 0;
    for(int j=0; j<M-1; ++j) sum_C_except_last += candC[j];
    
    for(int i = 0; i < N - 1; ++i) A[i][M - 1] = candR[i];
    for(int j = 0; j < M - 1; ++j) A[N - 1][j] = candC[j];
    A[N - 1][M - 1] = candR[N - 1] - sum_C_except_last;
    
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            cout << A[i][j] << (j == M - 1 ? "" : " ");
        }
        cout << "\n";
    }
    exit(0);
}

void solve_identical_rows() {
    if (N != M) {
        int K = abs(N - M);
        vector<vector<int>> prev_r(M + 1, vector<int>(K, -1));
        vector<vector<int>> choice(M + 1, vector<int>(K, -1));
        prev_r[0][0] = 0;
        
        for (int j = 0; j < M; ++j) {
            long long val = C[0][j] % K;
            for (int r = 0; r < K; ++r) {
                if (prev_r[j][r] != -1) {
                    int nxt1 = (r + val) % K;
                    prev_r[j + 1][nxt1] = r;
                    choice[j + 1][nxt1] = 1;
                    
                    int nxt2 = (r - val + K) % K;
                    prev_r[j + 1][nxt2] = r;
                    choice[j + 1][nxt2] = -1;
                }
            }
        }
        int curr = 0; 
        vector<int> s(M);
        for (int j = M; j > 0; --j) {
            s[j - 1] = choice[j][curr];
            curr = prev_r[j][curr];
        }
        
        long long S_C = 0;
        for (int j = 0; j < M; ++j) S_C += s[j] * C[0][j];
        long long X = S_C / (M - N);
        
        vector<long long> candR(N, X);
        vector<long long> candC(M);
        for (int j = 0; j < M; ++j) candC[j] = s[j] * C[0][j] + X;
        
        print_grid(candR, candC);
    } else {
        long long maxS = 0;
        for (int j = 0; j < M; ++j) maxS += C[0][j];
        if (maxS % 2 != 0) return; 
        long long target = maxS / 2;
        vector<int> dp(target + 1, -1);
        dp[0] = 0;
        for (int j = 0; j < M; ++j) {
            for (long long v = target; v >= C[0][j]; --v) {
                if (dp[v - C[0][j]] != -1 && dp[v] == -1) {
                    dp[v] = j + 1;
                }
            }
        }
        vector<int> s(M, -1);
        long long curr = target;
        while (curr > 0) {
            int j = dp[curr] - 1;
            s[j] = 1;
            curr -= C[0][j];
        }
        
        vector<long long> candR(N, 0);
        vector<long long> candC(M);
        for (int j = 0; j < M; ++j) candC[j] = s[j] * C[0][j];
        
        print_grid(candR, candC);
    }
}

void solve_identical_cols() {
    if (N != M) {
        int K = abs(N - M);
        vector<vector<int>> prev_r(N + 1, vector<int>(K, -1));
        vector<vector<int>> choice(N + 1, vector<int>(K, -1));
        prev_r[0][0] = 0;
        
        for (int i = 0; i < N; ++i) {
            long long val = C[i][0] % K;
            for (int r = 0; r < K; ++r) {
                if (prev_r[i][r] != -1) {
                    int nxt1 = (r + val) % K;
                    prev_r[i + 1][nxt1] = r;
                    choice[i + 1][nxt1] = 1;
                    
                    int nxt2 = (r - val + K) % K;
                    prev_r[i + 1][nxt2] = r;
                    choice[i + 1][nxt2] = -1;
                }
            }
        }
        int curr = 0;
        vector<int> s(N);
        for (int i = N; i > 0; --i) {
            s[i - 1] = choice[i][curr];
            curr = prev_r[i][curr];
        }
        
        long long S_R = 0;
        for (int i = 0; i < N; ++i) S_R += s[i] * C[i][0];
        long long X = -S_R / (N - M);
        
        vector<long long> candR(N);
        for (int i = 0; i < N; ++i) candR[i] = s[i] * C[i][0] + X;
        vector<long long> candC(M, X);
        
        print_grid(candR, candC);
    } else {
        long long maxS = 0;
        for (int i = 0; i < N; ++i) maxS += C[i][0];
        if (maxS % 2 != 0) return;
        long long target = maxS / 2;
        vector<int> dp(target + 1, -1);
        dp[0] = 0;
        for (int i = 0; i < N; ++i) {
            for (long long v = target; v >= C[i][0]; --v) {
                if (dp[v - C[i][0]] != -1 && dp[v] == -1) {
                    dp[v] = i + 1;
                }
            }
        }
        vector<int> s(N, -1);
        long long curr = target;
        while (curr > 0) {
            int i = dp[curr] - 1;
            s[i] = 1;
            curr -= C[i][0];
        }
        
        vector<long long> candR(N);
        for (int i = 0; i < N; ++i) candR[i] = s[i] * C[i][0];
        vector<long long> candC(M, 0);
        
        print_grid(candR, candC);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N >> M)) return 0;
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            cin >> C[i][j];
        }
    }
    
    if (N == 1 && M == 1) {
        cout << 0 << "\n";
        return 0;
    }

    bool all_zero = true;
    int i0 = -1, j0 = -1;
    long long D_i0_j0 = 0;
    
    for (int i = 1; i < N; ++i) {
        for (int j = 1; j < M; ++j) {
            long long d = C[i][0]*C[i][0] + C[0][j]*C[0][j] - C[i][j]*C[i][j] - C[0][0]*C[0][0];
            if (d != 0) {
                all_zero = false;
                i0 = i; j0 = j;
                D_i0_j0 = d;
                break;
            }
        }
        if (!all_zero) break;
    }

    if (all_zero) {
        bool identical_rows = true;
        for(int i=1; i<N; ++i) {
            for(int j=0; j<M; ++j) {
                if(C[i][j] != C[0][j]) identical_rows = false;
            }
        }
        if (identical_rows) {
            solve_identical_rows();
            return 0;
        }
        
        bool identical_cols = true;
        for(int i=0; i<N; ++i) {
            for(int j=1; j<M; ++j) {
                if(C[i][j] != C[i][0]) identical_cols = false;
            }
        }
        if (identical_cols) {
            solve_identical_cols();
            return 0;
        }
    } 
    else {
        long long c0_options[2] = { C[0][0], -C[0][0] };
        long long cj0_options[2] = { C[0][j0], -C[0][j0] };
        
        for (int opt_c0 = 0; opt_c0 < 2; ++opt_c0) {
            for (int opt_cj0 = 0; opt_cj0 < 2; ++opt_cj0) {
                long long c0 = c0_options[opt_c0];
                long long cj0 = cj0_options[opt_cj0];
                
                if (cj0 == c0) continue;
                if (D_i0_j0 % (2 * (cj0 - c0)) != 0) continue;
                
                long long r_i0 = D_i0_j0 / (2 * (cj0 - c0));
                if (r_i0 == 0) continue;
                
                bool possible = true;
                vector<long long> candR(N), candC(M);
                candC[0] = c0;
                
                for (int i = 0; i < N; ++i) {
                    long long num = C[i][0]*C[i][0] + C[0][j0]*C[0][j0] - C[i][j0]*C[i][j0] - C[0][0]*C[0][0];
                    if (num % (2 * (cj0 - c0)) != 0) { possible = false; break; }
                    candR[i] = num / (2 * (cj0 - c0));
                }
                if (!possible) continue;
                
                for (int j = 0; j < M; ++j) {
                    long long num = C[i0][0]*C[i0][0] + C[0][j]*C[0][j] - C[i0][j]*C[i0][j] - C[0][0]*C[0][0];
                    if (num % (2 * r_i0) != 0) { possible = false; break; }
                    candC[j] = c0 + num / (2 * r_i0);
                }
                if (!possible) continue;
                
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < M; ++j) {
                        if (abs(candR[i] - candC[j]) != C[i][j]) {
                            possible = false; break;
                        }
                    }
                    if (!possible) break;
                }
                if (!possible) continue;
                
                // 보정값 X 계산
                long long sumR = 0, sumC = 0;
                for (long long x : candR) sumR += x;
                for (long long x : candC) sumC += x;
                
                long long diff = sumC - sumR;
                if (N == M) {
                    if (diff != 0) possible = false;
                } else {
                    if (diff % (N - M) != 0) possible = false;
                }
                if (!possible) continue;
                
                long long X = 0;
                if (N != M) X = diff / (N - M);
                
                for (int i = 0; i < N; ++i) candR[i] += X;
                for (int j = 0; j < M; ++j) candC[j] += X;
                
                print_grid(candR, candC);
            }
        }
    }
    
    return 0;
}