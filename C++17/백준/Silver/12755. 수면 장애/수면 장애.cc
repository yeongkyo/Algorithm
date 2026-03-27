#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long N;
    cin >> N;

    long long len = 1;
    long long count = 9;
    long long start = 1;

    while (N > len * count) {
        N -= len * count;
        len++;
        count *= 10;
        start *= 10;
    }

    long long num = start + (N - 1) / len;
    string s = to_string(num);

    cout << s[(N - 1) % len] << '\n';

    return 0;
}