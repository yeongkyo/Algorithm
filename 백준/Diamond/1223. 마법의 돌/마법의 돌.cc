#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

long long dp[35][3][3][65][3];

long long solve(int L, int cL, int cR, int alts, int cmp, int n, int k, const vector<int>& S) {
    if (alts > k) return 0;
    
    int R = n - 1 - L;
    
    if (L > R) {
        if (cmp == 2) return 0;
        int extra = 0;
        if (n % 2 == 0 && cL != cR) extra = 1;
        if (alts + extra > k) return 0;
        return 1;
    }

    if (dp[L][cL][cR][alts][cmp] != -1) {
        return dp[L][cL][cR][alts][cmp];
    }

    long long sum = 0;
    int limit_L_start = 0, limit_L_end = 1;
    if (S[L] != -1) { 
        limit_L_start = limit_L_end = S[L]; 
    }
    
    int limit_R_start = 0, limit_R_end = 1;
    if (S[R] != -1) { 
        limit_R_start = limit_R_end = S[R]; 
    }

    for (int chL = limit_L_start; chL <= limit_L_end; ++chL) {
        for (int chR = limit_R_start; chR <= limit_R_end; ++chR) {
            if (L == R && chL != chR) continue;

            int ncmp = cmp;
            if (ncmp == 0) {
                if (chL < chR) ncmp = 1;
                else if (chL > chR) ncmp = 2;
            }

            int nalts = alts;
            if (cL != 2 && chL != cL) nalts++;
            if (cR != 2 && chR != cR) nalts++;

            sum += solve(L + 1, chL, chR, nalts, ncmp, n, k, S);
        }
    }

    return dp[L][cL][cR][alts][cmp] = sum;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, k;
    long long i;
    if (!(cin >> n >> k >> i)) return 0;

    vector<int> S(n, -1);

    memset(dp, -1, sizeof(dp));
    long long total = solve(0, 2, 2, 0, 0, n, k, S);

    if (total < i) {
        cout << "NO SUCH STONE\n";
        return 0;
    }

    for (int p = 0; p < n; ++p) {
        S[p] = 0;
        memset(dp, -1, sizeof(dp));
        long long cnt = solve(0, 2, 2, 0, 0, n, k, S);

        if (cnt >= i) {
            continue;
        } else {
            S[p] = 1;
            i -= cnt;
        }
    }

    for (int p = 0; p < n; ++p) {
        if (S[p] == 0) cout << 'I';
        else cout << 'X';
    }
    cout << "\n";

    return 0;
}