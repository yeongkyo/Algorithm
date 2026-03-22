#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int a, t;
    if (!(cin >> a >> t)) return 0;

    int p = 10 + 2 * (25 - a + t);

    if (p < 0) {
        cout << 0 << "\n";
    } else {
        cout << p << "\n";
    }

    return 0;
}