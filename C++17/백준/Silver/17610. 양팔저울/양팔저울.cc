#include <bits/stdc++.h>
using namespace std;

int k;
int w[13];
bool possible[2600001];

void dfs(int idx, int sum) {
    if (idx == k) {
        possible[abs(sum)] = true;
        return;
    }
    dfs(idx + 1, sum);
    dfs(idx + 1, sum + w[idx]);
    dfs(idx + 1, sum - w[idx]);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> k;
    int S = 0;
    for (int i = 0; i < k; i++) {
        cin >> w[i];
        S += w[i];
    }

    dfs(0, 0);

    int ans = 0;
    for (int i = 1; i <= S; i++) {
        if (!possible[i]) ans++;
    }

    cout << ans << '\n';
    return 0;
}