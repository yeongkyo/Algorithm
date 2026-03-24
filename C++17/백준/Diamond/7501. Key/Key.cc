#include <iostream>

using namespace std;

using u64 = unsigned long long;
using u128 = __uint128_t;

u64 power(u64 base, u64 exp, u64 mod) {
    u64 res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (u128)res * base % mod;
        base = (u128)base * base % mod;
        exp /= 2;
    }
    return res;
}

bool is_prime(u64 n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    u64 d = n - 1;
    int s = 0;
    while (d % 2 == 0) {
        d /= 2;
        s++;
    }

    static const u64 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 a : bases) {
        if (n <= a) break;
        u64 x = power(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; r++) {
            x = (u128)x * x % n;
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    u64 A, B;
    if (cin >> A >> B) {
        bool first = true;
        for (u64 K = A; K <= B; K++) {
            if (K % 2 == 0) continue;
            if (K == 9 || is_prime(K)) {
                if (!first) cout << " ";
                cout << K;
                first = false;
            }
        }
        cout << "\n";
    }
    return 0;
}