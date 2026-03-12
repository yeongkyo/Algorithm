#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int L;
    cin >> L;
    vector<int> S(L);
    for (int i = 0; i < L; i++) {
        cin >> S[i];
    }
    int n;
    cin >> n;

    sort(S.begin(), S.end());

    int left_bound = 0;
    int right_bound = 0;

    for (int s : S) {
        if (s == n) {
            cout << 0 << "\n";
            return 0;
        }
        if (s < n) {
            left_bound = s;
        } else {
            right_bound = s;
            break;
        }
    }

    int start = left_bound + 1;
    int end = right_bound - 1;

    int result = (end - n) + (n - start) * (end - n + 1);

    cout << result << "\n";

    return 0;
}