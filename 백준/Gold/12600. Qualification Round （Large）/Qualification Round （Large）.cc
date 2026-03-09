#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void solve(int tc) {
    int p;
    long long c;
    cin >> p >> c;
    
    vector<long long> s(p);
    for (int i = 0; i < p; ++i) {
        cin >> s[i];
    }

    long long low = 0;
    long long high = 6000000000000000000LL;
    long long ans = 0;

    while (low <= high) {
        long long mid = low + (high - low) / 2;
        __int128 total = 0;
        
        for (int i = 0; i < p; ++i) {
            total += min((long long)s[i], mid);
        }

        if (total >= (__int128)mid * c) {
            ans = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    
    cout << "Case #" << tc << ": " << ans << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int t;
    if (cin >> t) {
        for (int i = 1; i <= t; ++i) {
            solve(i);
        }
    }
    return 0;
}