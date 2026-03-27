#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N, M;
    cin >> N >> M;

    vector<vector<int>> A(N, vector<int>(M)), B(N, vector<int>(M));

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            cin >> A[i][j];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            cin >> B[i][j];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << A[i][j] + B[i][j];
            if (j + 1 < M) cout << ' ';
        }
        cout << '\n';
    }

    return 0;
}