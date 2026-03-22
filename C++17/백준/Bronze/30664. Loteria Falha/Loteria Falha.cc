#include <iostream>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s;
    while (cin >> s && s != "0") {
        int remainder = 0;
        for (char c : s) {
            remainder = (remainder * 10 + (c - '0')) % 42;
        }

        if (remainder == 0) {
            cout << "PREMIADO" << "\n";
        } else {
            cout << "TENTE NOVAMENTE" << "\n";
        }
    }

    return 0;
}