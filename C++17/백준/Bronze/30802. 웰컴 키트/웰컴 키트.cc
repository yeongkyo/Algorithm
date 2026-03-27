#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long N, T, P;
    cin >> N;

    long long a[6];
    for (int i = 0; i < 6; i++) cin >> a[i];

    cin >> T >> P;

    long long tshirt = 0;
    for (int i = 0; i < 6; i++) {
        tshirt += (a[i] + T - 1) / T;
    }

    cout << tshirt << '\n';
    cout << N / P << ' ' << N % P << '\n';

    return 0;
}