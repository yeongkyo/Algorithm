#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int N;
    if (!(cin >> N)) return 0;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; ++i) cin >> A[i];
    for (int i = 0; i < N; ++i) cin >> B[i];
    
    vector<pair<long long, long long>> v(N);
    for (int i = 0; i < N; ++i) v[i] = {B[i], A[i]};
    sort(v.begin(), v.end());
    
    long long prev = 0;
    long long ans = 0;
    int i = 0;
    while (i < N) {
        int j = i;
        while (j < N && v[j].first == v[i].first) ++j;
        long long group_max = prev;
        for (int k = i; k < j; ++k) {
            long long b = v[k].first;
            long long a = v[k].second;
            long long need = max(a, max(b, prev));
            long long diff = need - a;
            long long inc = diff > 0 ? (diff + 29) / 30 : 0;
            long long T = a + inc * 30;
            group_max = max(group_max, T);
            ans += inc;
        }
        prev = group_max;
        i = j;
    }
    
    cout << ans << '\n';
    return 0;
}