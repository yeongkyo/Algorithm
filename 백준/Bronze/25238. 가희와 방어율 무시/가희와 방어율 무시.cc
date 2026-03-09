#include <iostream>

using namespace std;

int main() {
    double a, b;
    if (!(cin >> a >> b)) return 0;

    double actual_armor = a * (100.0 - b) / 100.0;

    if (actual_armor < 100) {
        cout << 1 << endl;
    } else {
        cout << 0 << endl;
    }

    return 0;
}