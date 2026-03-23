#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    string first_word;
    cin >> first_word;

    vector<int> first_cnt(26, 0);
    for (char c : first_word) {
        first_cnt[c - 'A']++;
    }

    int ans = 0;
    for (int i = 1; i < n; i++) {
        string word;
        cin >> word;

        vector<int> cnt(26, 0);
        for (char c : word) {
            cnt[c - 'A']++;
        }

        int abs_diff = 0;
        for (int j = 0; j < 26; j++) {
            abs_diff += abs(first_cnt[j] - cnt[j]);
        }

        int len_diff = abs((int)first_word.length() - (int)word.length());

        if (len_diff == 0 && (abs_diff == 0 || abs_diff == 2)) {
            ans++;
        } else if (len_diff == 1 && abs_diff == 1) {
            ans++;
        }
    }

    cout << ans << "\n";

    return 0;
}