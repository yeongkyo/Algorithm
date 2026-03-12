#include <iostream>
#include <vector>

using namespace std;

int spf[100001];

void sieve(int n) {
    for (int i = 2; i <= n; ++i) spf[i] = i;
    for (int i = 2; i * i <= n; ++i) {
        if (spf[i] == i) {
            for (int j = i * i; j <= n; j += i)
                if (spf[j] == j) spf[j] = i;
        }
    }
}

long long count_p(long long n, long long p) {
    long long res = 0;
    while (n > 0) {
        res += n / p;
        n /= p;
    }
    return res;
}

bool check(long long n, long long k, int i) {
    int temp = i;
    while (temp > 1) {
        int p = spf[temp];
        int cnt = 0;
        while (temp % p == 0) {
            cnt++;
            temp /= p;
        }
        long long vp_comb = count_p(n, p) - count_p(k, p) - count_p(n - k, p);
        if (cnt > vp_comb) return false;
    }
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long S, F;
    int M;
    if (!(cin >> S >> F >> M)) return 0;

    sieve(M);

    long long n = S + F;
    long long k = S;

    for (int i = M; i >= 1; --i) {
        if (check(n, k, i)) {
            cout << i << "\n";
            return 0;
        }
    }

    cout << -1 << "\n";

    return 0;
}