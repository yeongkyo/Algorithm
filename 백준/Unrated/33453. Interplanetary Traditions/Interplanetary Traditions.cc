#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

const int MAXSQ = 225000;
long long val[MAXSQ * 2];
long long C_arr[MAXSQ * 2];
__int128_t S_p[MAXSQ * 2];
int ind1[MAXSQ];
int ind2[MAXSQ];

long long n;
long long sq;
vector<long long> primes;
vector<__int128_t> sum_primes;

int get_idx(long long v) {
    if (v <= sq) return ind1[v];
    else return ind2[n / v];
}

pair<__int128_t, long long> dfs(long long v, int j) {
    long long p_j = (j < primes.size()) ? primes[j] : sq + 1;
    if (v < p_j) return {0, 1};
    
    int idx = get_idx(v);
    __int128_t ans_res = S_p[idx] - sum_primes[j];
    long long F_res = 1 + C_arr[idx] - j;
    
    for (int k = j; k < primes.size(); ++k) {
        long long p_k = primes[k];
        if (p_k * p_k > v) break;
        
        auto child = dfs(v / p_k, k + 1);
        ans_res += (__int128_t)p_k * (child.second - 1);
        F_res += child.second - 1;
    }
    return {ans_res, F_res};
}

void print(__int128_t x) {
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
    
    if (!(cin >> n)) return 0;
    if (n == 1) {
        cout << 0 << "\n";
        return 0;
    }
    
    sq = sqrt(n);
    int tot = 0;
    
    for (long long l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l);
        long long v = n / l;
        val[++tot] = v;
        C_arr[tot] = v - 1;
        S_p[tot] = (__int128_t)v * (v + 1) / 2 - 1;
        if (v <= sq) ind1[v] = tot;
        else ind2[n / v] = tot;
    }
    
    sum_primes.push_back(0);
    vector<bool> is_prime(sq + 1, true);
    for (long long p = 2; p <= sq; ++p) {
        if (is_prime[p]) {
            primes.push_back(p);
            sum_primes.push_back(sum_primes.back() + p);
            for (long long i = p * p; i <= sq; i += p) {
                is_prime[i] = false;
            }
        }
    }
    
    for (size_t i = 0; i < primes.size(); ++i) {
        long long p = primes[i];
        long long p2 = p * p;
        __int128_t sp_prev = sum_primes[i];
        long long cp_prev = i;
        
        for (int j = 1; j <= tot; ++j) {
            long long v = val[j];
            if (v < p2) break;
            
            long long nv = v / p;
            int idx = get_idx(nv);
            
            C_arr[j] -= (C_arr[idx] - cp_prev);
            S_p[j] -= (__int128_t)p * (S_p[idx] - sp_prev);
        }
    }
    
    pair<__int128_t, long long> res = dfs(n, 0);
    
    __int128_t total_cost = (__int128_t)n - res.second + res.first;
    
    print(total_cost);
    
    return 0;
}