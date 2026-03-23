#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n;
    int k, q;
    cin >> n >> k >> q;

    while (q--) {
        long long x, y;
        cin >> x >> y;

        if (k == 1) {
            cout << abs(x - y) << "\n";
            continue;
        }

        long long dist = 0;
        while (x != y) {
            if (x > y) x = (x - 2) / k + 1;
            else y = (y - 2) / k + 1;
            dist++;
        }
        cout << dist << "\n";
    }

    return 0;
}