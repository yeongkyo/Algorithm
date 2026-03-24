#include <iostream>
#include <vector>

using namespace std;

const int MOD = 1000000007;
int n;
int tree[12][100005];

void add(int k, int idx, int val) {
    while (idx <= n) {
        tree[k][idx] = (tree[k][idx] + val) % MOD;
        idx += idx & -idx;
    }
}

int query(int k, int idx) {
    int sum = 0;
    while (idx > 0) {
        sum = (sum + tree[k][idx]) % MOD;
        idx -= idx & -idx;
    }
    return sum;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    for (int i = 0; i < n; i++) {
        int x;
        cin >> x;

        int counts[12] = {0};
        counts[1] = 1;

        for (int k = 2; k <= 11; k++) {
            counts[k] = query(k - 1, x - 1);
        }

        for (int k = 1; k <= 11; k++) {
            add(k, x, counts[k]);
        }
    }

    cout << query(11, n) << "\n";

    return 0;
}