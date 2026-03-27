#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long N;
    cin >> N;

    long long limit = 10;
    while (limit <= N) limit *= 10;

    int cnt = 0;
    while (N * 2 < limit) {
        N *= 2;
        cnt++;
    }

    cout << cnt;
    return 0;
}