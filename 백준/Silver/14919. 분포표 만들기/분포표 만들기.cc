#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int m;
    if (!(cin >> m)) return 0;

    vector<int> counts(m, 0);
    double val;
    while (cin >> val) {
        int index = (int)(val * m + 1e-9);
        if (index >= 0 && index < m) {
            counts[index]++;
        }
    }

    for (int i = 0; i < m; ++i) {
        if (i > 0) cout << " ";
        cout << counts[i];
    }
    cout << endl;

    return 0;
}