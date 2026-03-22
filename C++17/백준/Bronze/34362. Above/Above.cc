#include <iostream>
#include <iomanip>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    double n;
    if (!(cin >> n)) return 0;

    cout << fixed << setprecision(4) << n - 0.3 << "\n";

    return 0;
}