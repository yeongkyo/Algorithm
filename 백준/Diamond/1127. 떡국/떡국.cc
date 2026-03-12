#include<iostream>
#include<vector>
#include<cstring>
using namespace std;

int N, K, c[55][55], f[55][55], v[55], ans, S, T, m, x;
char A[55][55];
vector<int> g[55], o[55];

int D(int x) {
    if(x == T) return 1;
    v[x] = 1;
    for(int i = 0; i <= T; ++i)
        if(!v[i] && c[x][i] > f[x][i] && D(i)) return f[x][i]++, f[i][x]--, 1;
    return 0;
}

int main() {
    ios_base::sync_with_stdio(0); cin.tie(0);
    
    cin >> N; S = N; T = N + 1;
    for(int i = 0; i < N; ++i) cin >> A[i];
    cin >> K;
    
    for(int t = 0; t < 2; ++t) {
        for(int i = 0; i < N; ++i) {
            cin >> m;
            while(m--) cin >> x, (t ? o[i] : g[i]).push_back(x);
        }
    }
    
    for(int k = 0; k < K; ++k) {
        memset(c, 0, sizeof(c)); 
        memset(f, 0, sizeof(f));
        
        for(int i = 0; i < N; ++i) {
            int hg = 0, ho = 0;
            for(int j : g[i]) hg |= (j == k);
            for(int j : o[i]) ho |= (j == k);
            
            hg ? c[S][i] = 1e9 : (!ho ? c[i][T] = 1e9 : 0);
            
            for(int j = 0; j < N; ++j) if(A[i][j] == '1') c[i][j] = 1;
        }
        
        while(memset(v, 0, sizeof(v)), D(S)) ans++;
    }
    cout << ans << "\n";
}