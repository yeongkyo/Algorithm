#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int x;
    if (!(cin >> x)) return 0;

    if (x >= 6) {
        cout << "Success!" << "\n";
    } else {
        cout << "Oh My God!" << "\n";
    }

    return 0;
}