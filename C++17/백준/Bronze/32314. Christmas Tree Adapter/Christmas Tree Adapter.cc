#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int a, w, v;
    cin >> a >> w >> v;

    if (w / v >= a) {
        cout << 1 << "\n";
    } else {
        cout << 0 << "\n";
    }

    return 0;
}