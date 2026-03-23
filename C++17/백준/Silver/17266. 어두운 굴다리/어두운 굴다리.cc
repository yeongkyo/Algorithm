#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    cin >> n >> m;

    vector<int> x(m);
    for (int i = 0; i < m; i++) {
        cin >> x[i];
    }

    int max_h = x[0];

    for (int i = 1; i < m; i++) {
        int gap = x[i] - x[i - 1];
        max_h = max(max_h, (gap + 1) / 2);
    }

    max_h = max(max_h, n - x[m - 1]);

    cout << max_h << "\n";

    return 0;
}