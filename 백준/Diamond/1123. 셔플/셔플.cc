#include <iostream>
#include <vector>
#include <string>

using namespace std;

const long long MOD = 600921647;

typedef vector<vector<long long>> matrix;

matrix mul(const matrix &A, const matrix &B, int size) {
    matrix C(size, vector<long long>(size, 0));
    for(int i = 0; i < size; ++i) {
        for(int k = 0; k < size; ++k) {
            if(A[i][k] == 0) continue;
            for(int j = 0; j < size; ++j) {
                if(B[k][j] == 0) continue;
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % MOD;
            }
        }
    }
    return C;
}

matrix power(matrix A, long long p, int size) {
    matrix res(size, vector<long long>(size, 0));
    for(int i = 0; i < size; ++i) res[i][i] = 1;
    while(p > 0) {
        if(p & 1) res = mul(res, A, size);
        A = mul(A, A, size);
        p >>= 1;
    }
    return res;
}

long long cnt[15][15];
long long dp[15][15];
int N, M;
vector<string> T;

long long get_sum(long long X) {
    if (X < 0) return 0;
    if (X <= 9) {
        long long ans = 0;
        for(int t = 1; t <= X; ++t) {
            for(int g = 0; g < M; ++g) {
                ans = (ans + dp[t][g]) % MOD;
            }
        }
        return ans;
    }
    
    int size = 9 * M + 1;
    matrix W(size, vector<long long>(size, 0));
    for(int g = 0; g < M; ++g) {
        int r = g * 9 + 0;
        for(int gp = 0; gp < M; ++gp) {
            if(T[gp][g] == 'Y') {
                for(int L = 1; L <= 9; ++L) {
                    int c = gp * 9 + (L - 1);
                    W[r][c] = (W[r][c] + cnt[g][L]) % MOD;
                }
            }
        }
        for(int k = 1; k <= 8; ++k) {
            W[g * 9 + k][g * 9 + k - 1] = 1;
        }
    }
    
    int r_sum = 9 * M;
    W[r_sum][r_sum] = 1;
    for(int g = 0; g < M; ++g) {
        for(int gp = 0; gp < M; ++gp) {
            if(T[gp][g] == 'Y') {
                for(int L = 1; L <= 9; ++L) {
                    int c = gp * 9 + (L - 1);
                    W[r_sum][c] = (W[r_sum][c] + cnt[g][L]) % MOD;
                }
            }
        }
    }
    
    matrix W_p = power(W, X - 9, size);
    
    vector<long long> V(size, 0);
    for(int g = 0; g < M; ++g) {
        for(int k = 0; k <= 8; ++k) {
            V[g * 9 + k] = dp[9 - k][g];
        }
    }
    long long initial_sum = 0;
    for(int t = 1; t <= 9; ++t) {
        for(int g = 0; g < M; ++g) {
            initial_sum = (initial_sum + dp[t][g]) % MOD;
        }
    }
    V[r_sum] = initial_sum;
    
    long long final_sum = 0;
    for(int c = 0; c < size; ++c) {
        final_sum = (final_sum + W_p[r_sum][c] * V[c]) % MOD;
    }
    
    return final_sum;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N)) return 0;
    
    for(int i = 0; i < N; ++i) {
        int g, L;
        cin >> g >> L;
        cnt[g][L]++;
    }
    
    cin >> M;
    T.resize(M);
    for(int i = 0; i < M; ++i) {
        cin >> T[i];
    }
    
    long long A, B;
    cin >> A >> B;
    
    for(int t = 1; t <= 9; ++t) {
        for(int g = 0; g < M; ++g) {
            dp[t][g] = cnt[g][t];
            for(int L = 1; L <= 9; ++L) {
                if (t - L > 0) {
                    for(int gp = 0; gp < M; ++gp) {
                        if (T[gp][g] == 'Y') {
                            dp[t][g] = (dp[t][g] + cnt[g][L] * dp[t - L][gp]) % MOD;
                        }
                    }
                }
            }
        }
    }
    
    long long ans = (get_sum(B) - get_sum(A - 1) + MOD) % MOD;
    cout << ans << "\n";
    
    return 0;
}