#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;
    int n = (int)s.size();

    // Build Suffix Array (doubling)
    vector<int> sa(n), rnk(n), tmp(n);
    for (int i = 0; i < n; i++) {
        sa[i] = i;
        rnk[i] = s[i] - 'a';
    }

    for (int k = 1;; k <<= 1) {
        auto cmp = [&](int i, int j) {
            if (rnk[i] != rnk[j]) return rnk[i] < rnk[j];
            int ri = (i + k < n) ? rnk[i + k] : -1;
            int rj = (j + k < n) ? rnk[j + k] : -1;
            return ri < rj;
        };

        sort(sa.begin(), sa.end(), cmp);

        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        }
        rnk = tmp;
        if (rnk[sa[n - 1]] == n - 1) break;
    }

    // Build LCP array (Kasai)
    vector<int> inv(n);
    for (int i = 0; i < n; i++) inv[sa[i]] = i;

    long long sumLCP = 0;
    int h = 0;
    for (int i = 0; i < n; i++) {
        int pos = inv[i];
        if (pos == 0) continue;
        int j = sa[pos - 1];
        while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
        sumLCP += h;
        if (h > 0) h--;
    }

    long long total = 1LL * n * (n + 1) / 2;
    long long ans = total - sumLCP;
    cout << ans << "\n";
    return 0;
}
