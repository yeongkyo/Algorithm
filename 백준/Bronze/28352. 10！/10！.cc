#include <iostream>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    long long factorial = 1;
    for (int i = 1; i <= n; i++) {
        factorial *= i;
    }

    long long seconds_in_week = 7LL * 24 * 60 * 60;

    cout << factorial / seconds_in_week << "\n";

    return 0;
}