#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<string> s(n);
    for (int i = 0; i < n; i++) cin >> s[i];

    int k;
    cin >> k;

    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    int ans = 0;
    do {
        string t = "";
        for (int i = 0; i < n; i++) t += s[idx[i]];

        int len = t.length();
        vector<int> pi(len, 0);
        int j = 0;
        for (int i = 1; i < len; i++) {
            while (j > 0 && t[i] != t[j]) j = pi[j - 1];
            if (t[i] == t[j]) pi[i] = ++j;
        }

        int p = len - pi[len - 1];
        int matches = (len % p == 0) ? (len / p) : 1;

        if (matches == k) ans++;

    } while (next_permutation(idx.begin(), idx.end()));

    cout << ans << "\n";
    return 0;
}