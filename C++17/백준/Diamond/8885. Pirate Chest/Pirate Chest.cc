#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int a_val, b_val, m, n;
long long grid[505][505];
__int128 ans = 0;
long long S;

int L_arr[505], R_arr[505], st[505];
long long C[505];
long long max_D[505];

void solve(int max_r, int max_c) {
    for (int r1 = 0; r1 < m; ++r1) {
        for (int j = 0; j < n; ++j) C[j] = 2e18;
        
        int limit_r = min(m, r1 + max_r);
        for (int r2 = r1; r2 < limit_r; ++r2) {
            int h = r2 - r1 + 1;
            for (int j = 0; j < n; ++j) {
                if (grid[r2][j] < C[j]) C[j] = grid[r2][j];
            }
            
            int top = 0;
            for (int j = 0; j < n; ++j) {
                while (top > 0 && C[st[top - 1]] >= C[j]) top--;
                L_arr[j] = (top == 0) ? -1 : st[top - 1];
                st[top++] = j;
            }
            
            top = 0;
            for (int j = n - 1; j >= 0; --j) {
                while (top > 0 && C[st[top - 1]] >= C[j]) top--;
                R_arr[j] = (top == 0) ? n : st[top - 1];
                st[top++] = j;
            }
            
            for (int w = 1; w <= max_c; ++w) max_D[w] = 0;
            
            for (int j = 0; j < n; ++j) {
                int W = R_arr[j] - L_arr[j] - 1;
                if (W > max_c) W = max_c;
                if (C[j] > max_D[W]) max_D[W] = C[j];
            }
            
            for (int w = max_c - 1; w >= 1; --w) {
                if (max_D[w + 1] > max_D[w]) max_D[w] = max_D[w + 1];
            }
            
            for (int w = 1; w <= max_c; ++w) {
                long long D = max_D[w];
                if (D > 0) {
                    long long A = (long long)h * w;
                    __int128 vol = (__int128)A * ( ((__int128)D * S - 1) / (S - A) );
                    if (vol > ans) ans = vol;
                }
            }
        }
    }
}

void print(__int128 x) {
    if (x == 0) {
        cout << 0 << "\n";
        return;
    }
    string s;
    while (x > 0) {
        s += (char)('0' + (x % 10));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    cout << s << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> a_val >> b_val >> m >> n)) return 0;
    
    S = (long long)m * n;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> grid[i][j];
        }
    }
    
    solve(a_val, b_val);
    if (a_val != b_val) {
        solve(b_val, a_val);
    }
    
    print(ans);
    return 0;
}