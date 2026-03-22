#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int x;
    if (!(cin >> x)) return 0;

    if (x < 13) {
        cout << x << "\n";
    } else {
        cout << x + 1 << "\n";
    }

    return 0;
}