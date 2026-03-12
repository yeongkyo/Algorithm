#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <iostream>
#include <vector>

using namespace std;

alignas(32) int W[200008];
alignas(32) int x[200008];
alignas(32) long long ans[200008];

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
    
    vector<int> S;
    S.reserve(M);
    for (int j = 0; j < M; ++j) {
        int p;
        cin >> p;
        if (S.empty()) {
            S.push_back(p);
            continue;
        }
        if (S.back() == p) continue;
        
        while (S.size() >= 2) {
            long long a = S[S.size() - 2];
            long long b = S.back();
            if ((b - a > 0 && p - b >= 0) || (b - a < 0 && p - b <= 0)) {
                S.pop_back();
            } else {
                break;
            }
        }
        S.push_back(p);
    }

    int N_pad = (N + 15) & ~15;
    for (int i = N; i < N_pad; ++i) {
        W[i] = 0;
        x[i] = 0;
        ans[i] = 0;
    }

    int m_sz = S.size();
    for (int j = 0; j < m_sz; ++j) {
        int p = S[j];
        #pragma GCC ivdep
        for (int i = 0; i < N_pad; ++i) {
            int cx = x[i];
            int L = p - W[i];
            int diff1 = L > cx ? L - cx : 0;
            int diff2 = cx > p ? cx - p : 0;
            ans[i] += diff1 + diff2;
            x[i] = cx + diff1 - diff2;
        }
    }

    for (int i = 0; i < N; ++i) {
        cout << ans[i] << "\n";
    }
    
    return 0;
}