#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

typedef __int128_t int128;

int128 stringTo128(string s) {
    int128 res = 0;
    for (char c : s) res = res * 10 + (c - '0');
    return res;
}

string int128ToString(int128 n) {
    if (n == 0) return "0";
    string s = "";
    while (n > 0) {
        s += (char)((n % 10) + '0');
        n /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    string sa, sb, sr;
    if (!(cin >> sa >> sb >> sr)) return 0;

    int128 a = stringTo128(sa);
    int128 b = stringTo128(sb);
    int128 r = stringTo128(sr);

    if (a == 0 || b == 0) {
        cout << 0 << "\n";
        return 0;
    }

    if (b > r / a) {
        cout << "overflow" << "\n";
    } else {
        int128 result = a * b;
        cout << int128ToString(result) << "\n";
    }

    return 0;
}