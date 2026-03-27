#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    for (int i = 2, p; i <= N; i++) cin >> p;

    vector<long long> A(N);
    for (int i = 0; i < N; i++) cin >> A[i];

    sort(A.begin(), A.end(), greater<long long>());

    long long sum = 0;
    for (int i = 0; i < N; i++) {
        sum += A[i];
        cout << sum << '\n';
    }

    return 0;
}