#include <iostream>
#include <vector>
#include <random>

using namespace std;

const long long MOD = 1e9 + 7;

long long power(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

void solve() {
    int n, d;
    if (!(cin >> n >> d)) return;

    vector<vector<long long>> a(n + 1, vector<long long>(n + 1, 0));
    for (int i = 1; i <= n + 1; ++i) {
        for (int j = 0; j <= n + 1 - i; ++j) {
            cin >> a[i - 1][j];
        }
    }

    long long c = a[n][0];
    if (c == 0) {
        cout << "NO\n";
        return;
    }

    int q = -1;
    for (int i = 1; i <= n; ++i) {
        if (a[n - i][i] != 0) {
            q = i;
            break;
        }
    }
    if (q == -1) q = n;

    if (n % q != 0) {
        cout << "NO\n";
        return;
    }
    for (int i = 0; i <= n; ++i) {
        long long expected = (i % q == 0) ? c : 0;
        if (a[n - i][i] != expected) {
            cout << "NO\n";
            return;
        }
    }

    vector<long long> b(q + 1, 0);
    b[q] = 1;
    long long inv_c = power(c, MOD - 2);
    for (int k = q - 1; k >= 1; --k) {
        long long sum = 0;
        for (int j = k + 1; j <= q; ++j) {
            if (n - j >= 0) {
                sum = (sum + b[j] * a[n - j][k]) % MOD;
            }
        }
        b[k] = (sum * inv_c) % MOD;
    }

    int m = n / q;
    vector<vector<long long>> Qpow(m + 2, vector<long long>(n + q + 1, 0));
    Qpow[0][0] = 1;
    for (int k = 1; k <= m + 1; ++k) {
        for (int i = 0; i <= (k - 1) * q; ++i) {
            if (!Qpow[k - 1][i]) continue;
            for (int j = 1; j <= q; ++j) {
                if (!b[j]) continue;
                Qpow[k][i + j] = (Qpow[k][i + j] + Qpow[k - 1][i] * b[j]) % MOD;
            }
        }
    }

    vector<long long> V(n + q + 1, 0);
    for (int i = 0; i <= n; ++i) {
        if (!a[i][0]) continue;
        for (int j = 1; j <= q; ++j) {
            if (!b[j]) continue;
            V[i + j] = (V[i + j] + a[i][0] * b[j]) % MOD;
        }
    }

    vector<long long> f(m + 2, 0);
    for (int k = m + 1; k >= 1; --k) {
        f[k] = V[k * q];
        if (f[k] == 0) continue;
        for (int i = 0; i <= k * q; ++i) {
            V[i] = (V[i] - f[k] * Qpow[k][i]) % MOD;
            if (V[i] < 0) V[i] += MOD;
        }
    }

    vector<long long> R(n + q + 1, 0);
    for (int k = 1; k <= m + 1; ++k) {
        if (!f[k]) continue;
        for (int i = 0; i <= k * q; ++i) {
            R[i] = (R[i] + f[k] * Qpow[k][i]) % MOD;
        }
    }

    bool possible = true;
    mt19937_64 rng(1337);
    for (int iter = 0; iter < 15; ++iter) {
        long long X = rng() % (MOD - 1) + 1;
        long long Y = rng() % (MOD - 1) + 1;

        vector<long long> pX(n + q + 1, 1), pY(n + q + 1, 1);
        for (int i = 1; i <= n + q; ++i) {
            pX[i] = pX[i - 1] * X % MOD;
            pY[i] = pY[i - 1] * Y % MOD;
        }

        vector<long long> P_m(n + 1, 0);
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n - i; ++j) {
                if (!a[i][j]) continue;
                long long term = a[i][j] * pX[i] % MOD * pY[j] % MOD;
                P_m[i + j] = (P_m[i + j] + term) % MOD;
            }
        }

        vector<long long> Q_k(q + 1, 0);
        for (int k = 1; k <= q; ++k) {
            Q_k[k] = b[k] * (pX[k] - pY[k] + MOD) % MOD;
        }

        for (int D = d + 1; D <= n + q; ++D) {
            long long val = 0;
            for (int k = 1; k <= min(D, q); ++k) {
                int M = D - k;
                if (M <= n) {
                    val = (val + P_m[M] * Q_k[k]) % MOD;
                }
            }
            long long r_term = R[D] * (pX[D] - pY[D] + MOD) % MOD;
            val = (val - r_term + MOD) % MOD;
            if (val != 0) {
                possible = false;
                break;
            }
        }
        if (!possible) break;
    }

    if (!possible) {
        cout << "NO\n";
    } else {
        cout << "YES\n";
        int r = n + q;
        cout << q << " " << r << "\n";
        for (int i = 0; i <= q; ++i) {
            cout << b[i] << (i == q ? "" : " ");
        }
        cout << "\n";
        for (int i = 0; i <= r; ++i) {
            cout << R[i] << (i == r ? "" : " ");
        }
        cout << "\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int t;
    if (cin >> t) {
        while (t--) solve();
    }
    return 0;
}