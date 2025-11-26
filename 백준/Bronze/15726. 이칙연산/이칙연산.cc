#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long A, B, C;
    if (cin >> A >> B >> C) {
        long long val1 = A * B / C;
        long long val2 = A * C / B;

        cout << max(val1, val2);
    }

    return 0;
}