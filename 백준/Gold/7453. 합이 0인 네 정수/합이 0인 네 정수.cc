#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;

    vector<int> A(n), B(n), C(n), D(n);
    for (int i = 0; i < n; i++) {
        cin >> A[i] >> B[i] >> C[i] >> D[i];
    }

    vector<int> AB;
    AB.reserve(1LL * n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AB.push_back(A[i] + B[j]);
        }
    }

    sort(AB.begin(), AB.end());

    long long ans = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int target = -(C[i] + D[j]);
            auto range = equal_range(AB.begin(), AB.end(), target);
            ans += (long long)(range.second - range.first);
        }
    }

    cout << ans << "\n";
    return 0;
}