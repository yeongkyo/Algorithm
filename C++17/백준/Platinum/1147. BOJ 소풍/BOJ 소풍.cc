#include <iostream>
#include <algorithm>

using namespace std;

long long gcd(long long a, long long b) {
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n, m;
    if (!(cin >> n >> m)) return 0;

    long long count = 0;

    for (long long X1 = 1; X1 <= n; ++X1) {
        for (long long Y1 = 0; Y1 <= m; ++Y1) {
            if (Y1 == 0) {
                long long X2 = 0;
                for (long long Y2 = 1; Y2 <= m; ++Y2) {
                    if (X1 % 2 == X2 % 2 && Y1 % 2 == Y2 % 2) {
                        long long W = max(X1, X2);
                        long long H = max(Y1, Y2);
                        if (W <= n && H <= m) {
                            count += (n - W + 1) * (m - H + 1);
                        }
                    }
                }
            } else {
                long long g = gcd(X1, Y1);
                long long u = X1 / g;
                long long v = Y1 / g;

                long long max_k = min(n / v, m / u);
                for (long long k = 1; k <= max_k; ++k) {
                    long long X2 = k * v;
                    long long Y2 = k * u;

                    if (X1 % 2 == X2 % 2 && Y1 % 2 == Y2 % 2) {
                        long long W = max(X1, X2);
                        long long H = max(Y1, Y2);
                        if (W <= n && H <= m) {
                            count += (n - W + 1) * (m - H + 1);
                        }
                    }
                }
            }
        }
    }

    cout << count << "\n";
    return 0;
}