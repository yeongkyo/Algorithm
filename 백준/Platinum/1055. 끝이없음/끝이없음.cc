#include <iostream>
#include <string>
#include <vector>

using namespace std;

string A, S;
long long N, min_val, max_val;
vector<long long> L;

char solve(int n, long long p) {
    if (n == 0) return A[p - 1];
    for (char c : S) {
        if (c == '$') {
            if (p <= L[n - 1]) {
                return solve(n - 1, p);
            } else {
                p -= L[n - 1];
            }
        } else {
            if (p == 1) {
                return c;
            } else {
                p--;
            }
        }
    }
    return '-';
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> A >> S >> N >> min_val >> max_val)) return 0;

    long long C = 0, K = 0;
    for (char c : S) {
        if (c == '$') C++;
        else K++;
    }

    if (C == 1) {
        long long total_len = A.length() + N * K;
        string text = "";
        for (int i = 1; i < S.length(); i++) {
            text += S[i];
        }

        for (long long p = min_val; p <= max_val; p++) {
            if (p > total_len) {
                cout << '-';
            } else if (p <= A.length()) {
                cout << A[p - 1];
            } else {
                long long idx = (p - A.length() - 1) % K;
                cout << text[idx];
            }
        }
        cout << '\n';
        return 0;
    }

    L.push_back(A.length());
    long long limit = N;
    for (int i = 1; i <= N; i++) {
        long long nxt = L.back() * C + K;
        L.push_back(nxt);
        if (nxt > 2000000000LL) {
            limit = i;
            break;
        }
    }

    for (long long p = min_val; p <= max_val; p++) {
        if (p > L[limit]) {
            cout << '-';
        } else {
            cout << solve(limit, p);
        }
    }
    cout << '\n';

    return 0;
}