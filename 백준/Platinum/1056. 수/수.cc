#include <iostream>
#include <unordered_map>

using namespace std;

unordered_map<long long, long long> memo;

long long get_v(long long x, int k) {
    long long l = 1, r = 1000000000LL; 
    long long ans = 1;
    while (l <= r) {
        long long mid = l + (r - l) / 2;
        long long temp = 1;
        bool ok = true;
        for (int i = 0; i < k; ++i) {
            if (temp > x / mid) {
                ok = false;
                break;
            }
            temp *= mid;
        }
        if (ok) {
            ans = mid;
            l = mid + 1;
        } else {
            r = mid - 1;
        }
    }
    return ans;
}

long long solve(long long x) {
    if (x == 1) return 0;
    if (x == 2) return 1;
    if (memo.count(x)) return memo[x];
    
    long long res = x - 1; 
    
    for (int k = 2; k <= 60; ++k) {
        long long v = get_v(x, k);
        
        if (v > 1) {
            long long p1 = 1;
            for (int i = 0; i < k; ++i) p1 *= v;
            long long cost1 = solve(v) + (x - p1) + 1;
            if (res > cost1) res = cost1;
        }
        
        long long p2 = 1;
        bool overflow = false;
        for (int i = 0; i < k; ++i) {
            if (p2 > (x + res) / (v + 1)) {
                overflow = true;
                break;
            }
            p2 *= (v + 1);
        }
        
        if (!overflow) {
            long long cost2 = solve(v + 1) + (p2 - x) + 1;
            if (res > cost2) res = cost2;
        }
        
        if (v == 1 && overflow) {
            break;
        }
    }
    
    return memo[x] = res;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    long long N;
    if (cin >> N) {
        cout << solve(N) << '\n';
    }
    
    return 0;
}