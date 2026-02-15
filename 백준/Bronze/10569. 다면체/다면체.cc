#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;
    while (T--) {
        int V, E;
        cin >> V >> E;
        cout << (E - V + 2) << "\n";
    }
    return 0;
}
