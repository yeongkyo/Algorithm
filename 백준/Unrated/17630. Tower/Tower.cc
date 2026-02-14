#include <bits/stdc++.h>
using namespace std;

static inline bool isPowerOfTwo(unsigned long long x) {
    return x && ((x & (x - 1ULL)) == 0ULL);
}

// ceil(log2(q)) for q>=2
static inline int ceilLog2(unsigned long long q) {
    // ceil(log2(q)) == bit_length(q-1)
    unsigned long long x = q - 1ULL;
    return 64 - __builtin_clzll(x);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;
    while (T--) {
        unsigned long long q;
        cin >> q;

        if (q == 1ULL) {
            cout << 0 << "\n";
            continue;
        }

        int k = ceilLog2(q);                 // minimal steps will be k+1
        bool pow2 = isPowerOfTwo(q);
        unsigned long long N = pow2 ? q : (q + 1ULL);

        // Build totals S[0..k] where S[k]=N and S[i]=ceil(S[i+1]/2)
        vector<unsigned long long> S(k + 1);
        S[k] = N;
        for (int i = k - 1; i >= 0; --i) {
            S[i] = (S[i + 1] + 1ULL) / 2ULL; // ceil half
        }
        // S[0] should be 1
        // Steps to build total N in k moves:
        // at move t (1..k), current length = t, use (1,t) if S[t] even else (2,t).
        vector<pair<int,int>> ops;
        ops.reserve(k + 1);

        for (int t = 1; t <= k; ++t) {
            int u = t;
            int l = (S[t] % 2ULL == 0ULL) ? 1 : 2;
            ops.push_back({l, u});
        }

        // Final move to write q on top:
        // if pow2: sum[1..k+1] = q (because total is q)
        // else: sum[2..k+1] = (q+1) - 1 = q
        ops.push_back({pow2 ? 1 : 2, k + 1});

        cout << ops.size() << "\n";
        for (auto [l, u] : ops) {
            cout << l << " " << u << "\n";
        }
    }
    return 0;
}
