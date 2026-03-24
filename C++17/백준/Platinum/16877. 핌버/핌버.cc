#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> P(N);
    int max_p = 0;
    for (int i = 0; i < N; i++) {
        cin >> P[i];
        if (P[i] > max_p) max_p = P[i];
    }

    vector<int> fib;
    int a = 1, b = 2;
    fib.push_back(1);
    fib.push_back(2);
    while (true) {
        int next_fib = a + b;
        if (next_fib > 3000000) break;
        fib.push_back(next_fib);
        a = b;
        b = next_fib;
    }

    vector<int> G(max_p + 1, 0);
    for (int i = 1; i <= max_p; i++) {
        unsigned long long check = 0;
        for (int f : fib) {
            if (i >= f) {
                check |= (1ULL << G[i - f]);
            } else {
                break;
            }
        }
        G[i] = __builtin_ctzll(~check);
    }

    int result = 0;
    for (int i = 0; i < N; i++) {
        result ^= G[P[i]];
    }

    if (result > 0) {
        cout << "koosaga\n";
    } else {
        cout << "cubelover\n";
    }

    return 0;
}