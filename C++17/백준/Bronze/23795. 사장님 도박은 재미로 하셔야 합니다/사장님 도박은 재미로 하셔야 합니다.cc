#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long total_lost = 0;
    int money;

    while (cin >> money && money != -1) {
        total_lost += money;
    }

    cout << total_lost << "\n";

    return 0;
}