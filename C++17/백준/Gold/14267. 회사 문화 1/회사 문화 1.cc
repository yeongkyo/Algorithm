#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    if (!(cin >> n >> m)) return 0;

    vector<int> parent(n + 1);
    vector<long long> score(n + 1, 0);

    for (int i = 1; i <= n; ++i) {
        cin >> parent[i];
    }

    for (int k = 0; k < m; ++k) {
        int i, w;
        cin >> i >> w;
        score[i] += w;
    }

    for (int i = 2; i <= n; ++i) {
        if (parent[i] != -1) {
            score[i] += score[parent[i]];
        }
    }

    for (int i = 1; i <= n; ++i) {
        cout << score[i] << (i == n ? "" : " ");
    }
    cout << "\n";

    return 0;
}