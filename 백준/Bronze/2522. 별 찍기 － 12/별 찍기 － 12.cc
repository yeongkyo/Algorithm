#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    // upper part: 1..N
    for (int i = 1; i <= N; i++) {
        cout << string(N - i, ' ') << string(i, '*') << "\n";
    }
    // lower part: N-1..1
    for (int i = N - 1; i >= 1; i--) {
        cout << string(N - i, ' ') << string(i, '*') << "\n";
    }
    return 0;
}
