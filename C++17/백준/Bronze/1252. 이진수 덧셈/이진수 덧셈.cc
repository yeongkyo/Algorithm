#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string a, b;
    cin >> a >> b;

    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());

    string res;
    int carry = 0;
    int n = max(a.size(), b.size());

    for (int i = 0; i < n; i++) {
        int x = (i < (int)a.size() ? a[i] - '0' : 0);
        int y = (i < (int)b.size() ? b[i] - '0' : 0);

        int sum = x + y + carry;
        res.push_back((sum % 2) + '0');
        carry = sum / 2;
    }

    if (carry) res.push_back('1');

    while (res.size() > 1 && res.back() == '0') res.pop_back();

    reverse(res.begin(), res.end());
    cout << res << '\n';

    return 0;
}