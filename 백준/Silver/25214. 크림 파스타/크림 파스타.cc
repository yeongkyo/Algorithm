#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    int min_val;
    int max_diff = 0;

    for (int i = 0; i < n; i++) {
        int current;
        cin >> current;

        if (i == 0) {
            min_val = current;
            cout << 0 << " ";
        } else {
            max_diff = max(max_diff, current - min_val);
            
            min_val = min(min_val, current);
            
            cout << max_diff << " ";
        }
    }
    cout << "\n";

    return 0;
}