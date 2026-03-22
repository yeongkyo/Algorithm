#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n, c, p;
    if (!(cin >> n >> c >> p)) return 0;

    if (n <= c * p) {
        cout << "yes" << "\n";
    } else {
        cout << "no" << "\n";
    }

    return 0;
}