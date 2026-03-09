#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int K;
    cin >> K;
    string s;
    cin >> s;

    vector<int> counts;
    int cnt = 1;
    for (int i = 1; i < K; ++i) {
        if (s[i] == s[i - 1]) {
            cnt++;
        } else {
            counts.push_back(cnt);
            cnt = 1;
        }
    }
    counts.push_back(cnt);

    int maxLen = 0;
    if (counts.size() > 1) {
        for (int i = 0; i < (int)counts.size() - 1; ++i) {
            int currentMagnetLen = 2 * min(counts[i], counts[i + 1]);
            if (currentMagnetLen > maxLen) {
                maxLen = currentMagnetLen;
            }
        }
    }

    cout << maxLen << "\n";

    return 0;
}