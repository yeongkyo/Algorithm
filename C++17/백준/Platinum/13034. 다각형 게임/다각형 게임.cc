#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> G(N + 1, 0);

    for (int i = 2; i <= N; i++) {
        vector<bool> check(N + 1, false);
        for (int j = 0; j <= i - 2; j++) {
            int val = G[j] ^ G[i - 2 - j];
            if (val <= N) {
                check[val] = true;
            }
        }
        
        int mex = 0;
        while (check[mex]) {
            mex++;
        }
        G[i] = mex;
    }

    if (G[N] > 0) {
        cout << 1 << "\n";
    } else {
        cout << 2 << "\n";
    }

    return 0;
}