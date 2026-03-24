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
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    reverse(a.begin(), a.end());

    vector<int> pi(n, 0);
    int j = 0;
    int max_pi = 0;
    int best_i = 0;

    for (int i = 1; i < n; i++) {
        while (j > 0 && a[i] != a[j]) {
            j = pi[j - 1];
        }
        if (a[i] == a[j]) {
            j++;
            pi[i] = j;
        }
    }

    for (int i = 0; i < n; i++) {
        if (pi[i] > max_pi) {
            max_pi = pi[i];
            best_i = i;
        }
    }

    int L = best_i + 1;
    int k = n - L;
    int p = L - max_pi;

    cout << k << " " << p << "\n";

    return 0;
}