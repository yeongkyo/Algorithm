#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    int sum = 0, cur = 0;

    for (char c : s) {
        if (c == ',') {
            sum += cur;
            cur = 0;
        } else {
            cur = cur * 10 + (c - '0');
        }
    }

    sum += cur;

    cout << sum << '\n';
    return 0;
}