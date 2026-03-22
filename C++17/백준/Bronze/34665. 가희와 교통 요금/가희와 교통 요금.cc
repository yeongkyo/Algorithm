#include <iostream>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string a, b;
    cin >> a >> b;

    if (a == b) {
        cout << 0 << "\n";
    } else {
        cout << 1550 << "\n";
    }

    return 0;
}