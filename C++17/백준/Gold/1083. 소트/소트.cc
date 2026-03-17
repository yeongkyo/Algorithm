#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> a(n);
    for (int i = 0; i < n; ++i) {
        cin >> a[i];
    }

    int s;
    cin >> s;

    for (int i = 0; i < n && s > 0; ++i) {
        int max_idx = i;
        for (int j = i + 1; j < n && j <= i + s; ++j) {
            if (a[j] > a[max_idx]) {
                max_idx = j;
            }
        }

        for (int j = max_idx; j > i; --j) {
            swap(a[j], a[j - 1]);
        }
        s -= (max_idx - i);
    }

    for (int i = 0; i < n; ++i) {
        cout << a[i] << (i == n - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}