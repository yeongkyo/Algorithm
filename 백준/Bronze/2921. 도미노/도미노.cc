#include <iostream>

using namespace std;

int main() {
    int n;
    cin >> n;

    long long total = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = i; j <= n; ++j) {
            total += (i + j);
        }
    }

    cout << total << endl;

    return 0;
}