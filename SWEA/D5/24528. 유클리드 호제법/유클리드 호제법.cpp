#include <iostream>

using namespace std;

void solve() {
    long long X, Y;
    cin >> X >> Y;

    long long p = Y;
    while (p <= X) {
        p *= Y;
    }

    long long b = (p + X - 1) / X;
    long long B = b * X;
    long long A = B * Y + X;

    cout << A << " " << B << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        while (T--) {
            solve();
        }
    }
    return 0;
}