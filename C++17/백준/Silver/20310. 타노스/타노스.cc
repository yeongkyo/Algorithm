#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s;
    cin >> s;

    int cnt0 = 0, cnt1 = 0;
    for (char c : s) {
        if (c == '0') cnt0++;
        else cnt1++;
    }

    int rem0 = cnt0 / 2;
    int rem1 = cnt1 / 2;

    vector<bool> deleted(s.length(), false);

    int deleted1 = 0;
    for (int i = 0; i < s.length() && deleted1 < rem1; i++) {
        if (s[i] == '1') {
            deleted[i] = true;
            deleted1++;
        }
    }

    int deleted0 = 0;
    for (int i = (int)s.length() - 1; i >= 0 && deleted0 < rem0; i--) {
        if (s[i] == '0') {
            deleted[i] = true;
            deleted0++;
        }
    }

    for (int i = 0; i < s.length(); i++) {
        if (!deleted[i]) cout << s[i];
    }
    cout << "\n";

    return 0;
}