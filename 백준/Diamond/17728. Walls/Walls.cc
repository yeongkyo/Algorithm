#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

alignas(32) int W[200008];
alignas(32) int x[200008];
alignas(32) long long ans[200008];
int P[200008];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int N, M;
    if (!(cin >> N >> M)) return 0;
    
    for (int i = 0; i < N; ++i) {
        int a, b;
        cin >> a >> b;
        W[i] = b - a;
        x[i] = a;
        ans[i] = 0;
    }
    
    for (int j = 0; j < M; ++j) {
        cin >> P[j];
    }

    int N_pad = (N + 7) & ~7;
    for(int i = N; i < N_pad; ++i) {
        W[i] = 0;
        x[i] = 0;
        ans[i] = 0;
    }

    const int BLOCK_SIZE = 1024;
    for (int i_start = 0; i_start < N_pad; i_start += BLOCK_SIZE) {
        int i_end = min(N_pad, i_start + BLOCK_SIZE);
        for (int j = 0; j < M; ++j) {
            int p = P[j];
            #pragma GCC ivdep
            #pragma GCC unroll 8
            for (int i = i_start; i < i_end; ++i) {
                int cx = x[i];
                int L = p - W[i];
                int diff1 = max(0, L - cx);
                int diff2 = max(0, cx - p);
                ans[i] += diff1 + diff2;
                x[i] = cx + diff1 - diff2;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        cout << ans[i] << "\n";
    }
    
    return 0;
}