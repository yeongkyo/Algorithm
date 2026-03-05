#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int K;
    if (!(cin >> K)) return 0;

    while (K--) {
        long long N, M;
        cin >> N >> M;

        long long min_val = min(N, M);
        long long result = 2 * (min_val - 1);

        cout << result << "\n";
    }

    return 0;
}