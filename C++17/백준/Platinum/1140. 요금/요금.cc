#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long T, K1, P1, K2, P2;
    if (!(cin >> T >> K1 >> P1 >> K2 >> P2)) return 0;

    P1 = min(P1, K1 * 10);
    P2 = min(P2, K2 * 10);

    long long ans = -1;
    auto update = [&](long long val) {
        if (ans == -1 || val < ans) ans = val;
    };

    long long limit = 2000000;
    long long max_x = T / K1 + 1;
    long long max_y = T / K2 + 1;

    for (long long x = 0; x <= min(limit, max_x); ++x) {
        long long rem = T - x * K1;
        long long c = x * P1;
        if (rem <= 0) {
            update(c);
        } else {
            long long k = rem / K2;
            long long rem2 = rem % K2;
            long long extra = min(rem2 * 10, P2);
            update(c + k * P2 + extra);
        }
    }

    for (long long x = max(0LL, max_x - limit); x <= max_x; ++x) {
        long long rem = T - x * K1;
        long long c = x * P1;
        if (rem <= 0) {
            update(c);
        } else {
            long long k = rem / K2;
            long long rem2 = rem % K2;
            long long extra = min(rem2 * 10, P2);
            update(c + k * P2 + extra);
        }
    }

    for (long long y = 0; y <= min(limit, max_y); ++y) {
        long long rem = T - y * K2;
        long long c = y * P2;
        if (rem <= 0) {
            update(c);
        } else {
            long long k = rem / K1;
            long long rem2 = rem % K1;
            long long extra = min(rem2 * 10, P1);
            update(c + k * P1 + extra);
        }
    }

    for (long long y = max(0LL, max_y - limit); y <= max_y; ++y) {
        long long rem = T - y * K2;
        long long c = y * P2;
        if (rem <= 0) {
            update(c);
        } else {
            long long k = rem / K1;
            long long rem2 = rem % K1;
            long long extra = min(rem2 * 10, P1);
            update(c + k * P1 + extra);
        }
    }

    cout << ans << "\n";
    return 0;
}