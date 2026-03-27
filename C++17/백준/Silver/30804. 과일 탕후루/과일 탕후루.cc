#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    vector<int> a(N);
    for (int i = 0; i < N; i++) cin >> a[i];

    int cnt[10] = {0};
    int kind = 0, left = 0, ans = 0;

    for (int right = 0; right < N; right++) {
        if (cnt[a[right]] == 0) kind++;
        cnt[a[right]]++;

        while (kind > 2) {
            cnt[a[left]]--;
            if (cnt[a[left]] == 0) kind--;
            left++;
        }

        ans = max(ans, right - left + 1);
    }

    cout << ans << '\n';
    return 0;
}