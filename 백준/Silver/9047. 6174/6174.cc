#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

int get_next(int n) {
    vector<int> digits(4);
    for (int i = 0; i < 4; ++i) {
        digits[i] = n % 10;
        n /= 10;
    }

    sort(digits.begin(), digits.end());
    int min_val = 0;
    for (int i = 0; i < 4; ++i) {
        min_val = min_val * 10 + digits[i];
    }

    sort(digits.begin(), digits.end(), greater<int>());
    int max_val = 0;
    for (int i = 0; i < 4; ++i) {
        max_val = max_val * 10 + digits[i];
    }

    return max_val - min_val;
}

void solve() {
    int n;
    cin >> n;

    int count = 0;
    while (n != 6174) {
        n = get_next(n);
        count++;
    }

    cout << count << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t;
    if (cin >> t) {
        while (t--) {
            solve();
        }
    }

    return 0;
}