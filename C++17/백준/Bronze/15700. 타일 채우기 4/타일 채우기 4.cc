#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n, m;
    if (!(cin >> n >> m)) return 0;

    cout << (n * m) / 2 << "\n";

    return 0;
}