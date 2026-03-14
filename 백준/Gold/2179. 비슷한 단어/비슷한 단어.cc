#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int getLCP(const string& a, const string& b) {
    int len = min(a.length(), b.length());
    for (int i = 0; i < len; i++) {
        if (a[i] != b[i]) return i;
    }
    return len;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    vector<string> orig(n);
    vector<pair<string, int>> sorted_words(n);

    for (int i = 0; i < n; i++) {
        cin >> orig[i];
        sorted_words[i] = {orig[i], i};
    }

    sort(sorted_words.begin(), sorted_words.end());

    int maxLCP = -1;
    for (int i = 0; i < n - 1; i++) {
        int lcp = getLCP(sorted_words[i].first, sorted_words[i + 1].first);
        if (lcp > maxLCP) {
            maxLCP = lcp;
        }
    }

    for (int i = 0; i < n; i++) {
        int sorted_idx = lower_bound(sorted_words.begin(), sorted_words.end(), make_pair(orig[i], i)) - sorted_words.begin();
        
        bool can_be_S = false;
        if (sorted_idx > 0 && getLCP(orig[i], sorted_words[sorted_idx - 1].first) >= maxLCP) can_be_S = true;
        if (sorted_idx < n - 1 && getLCP(orig[i], sorted_words[sorted_idx + 1].first) >= maxLCP) can_be_S = true;

        if (can_be_S) {
            for (int j = i + 1; j < n; j++) {
                if (getLCP(orig[i], orig[j]) == maxLCP) {
                    cout << orig[i] << "\n" << orig[j] << "\n";
                    return 0;
                }
            }
        }
    }

    return 0;
}