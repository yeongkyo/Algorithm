#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<int> A(10), B(10);
    for (int i = 0; i < 10; i++) cin >> A[i];
    for (int i = 0; i < 10; i++) cin >> B[i];

    int scoreA = 0, scoreB = 0;
    char last = 'D';

    for (int i = 0; i < 10; i++) {
        if (A[i] > B[i]) {
            scoreA += 3;
            last = 'A';
        } else if (A[i] < B[i]) {
            scoreB += 3;
            last = 'B';
        } else {
            scoreA += 1;
            scoreB += 1;
        }
    }

    cout << scoreA << ' ' << scoreB << '\n';

    if (scoreA > scoreB) cout << 'A';
    else if (scoreA < scoreB) cout << 'B';
    else cout << last;

    cout << '\n';
    return 0;
}