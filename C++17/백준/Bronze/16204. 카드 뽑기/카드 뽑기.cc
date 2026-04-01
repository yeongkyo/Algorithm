#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M, K;
    cin >> N >> M >> K;

    int same_O = min(M, K);
    int same_X = min(N - M, N - K);

    cout << same_O + same_X << "\n";

    return 0;
}