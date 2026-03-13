#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

bool isLucky(const string& s) {
    for (int i = 0; i < s.length() - 1; i++) {
        if (s[i] == s[i + 1]) return false;
    }
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s;
    cin >> s;

    sort(s.begin(), s.end());

    int luckyCount = 0;

    do {
        if (isLucky(s)) {
            luckyCount++;
        }
    } while (next_permutation(s.begin(), s.end()));

    cout << luckyCount << "\n";

    return 0;
}