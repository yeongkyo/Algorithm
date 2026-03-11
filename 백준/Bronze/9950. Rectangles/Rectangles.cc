#include <iostream>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int l, w, a;

    while (cin >> l >> w >> a && (l != 0 || w != 0 || a != 0)) {
        if (l == 0) {
            l = a / w;
        } else if (w == 0) {
            w = a / l;
        } else if (a == 0) {
            a = l * w;
        }

        cout << l << " " << w << " " << a << "\n";
    }

    return 0;
}