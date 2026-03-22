#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long m, k, c;
    if (!(cin >> m >> k >> c)) return 0;

    cout << m * c << "\n";

    return 0;
}