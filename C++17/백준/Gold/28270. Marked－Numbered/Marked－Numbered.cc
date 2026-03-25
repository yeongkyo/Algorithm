#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int c;
    if (!(cin >> c)) return 0;

    vector<int> m(c);
    for (int i = 0; i < c; i++) {
        cin >> m[i];
    }

    vector<int> n(c);
    vector<int> cnt(1000005, 0); 
    int prev_depth = 0;

    for (int i = 0; i < c; i++) {
        int current_depth = m[i];

        if (current_depth > prev_depth + 1) {
            cout << -1 << "\n";
            return 0;
        }

        cnt[current_depth]++;

        cnt[current_depth + 1] = 0;

        n[i] = cnt[current_depth];
        
        prev_depth = current_depth;
    }

    for (int i = 0; i < c; i++) {
        cout << n[i] << (i == c - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}