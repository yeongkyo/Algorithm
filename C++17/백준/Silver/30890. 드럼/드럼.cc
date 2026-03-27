#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int X, Y;
    cin >> X >> Y;

    int i = 1, j = 1;

    while (i <= X || j <= Y) {
        if (i > X) {
            cout << 2;
            j++;
        } else if (j > Y) {
            cout << 1;
            i++;
        } else {
            long long left = (long long)i * Y;
            long long right = (long long)j * X;

            if (left == right) {
                cout << 3;
                i++;
                j++;
            } else if (left < right) {
                cout << 1;
                i++;
            } else {
                cout << 2;
                j++;
            }
        }
    }

    return 0;
}