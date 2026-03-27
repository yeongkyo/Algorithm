#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;

    while (T--) {
        int N;
        cin >> N;

        vector<int> A(N + 1);
        for (int i = 1; i <= N; i++) cin >> A[i];

        vector<bool> visited(N + 1, false);
        int cur = 1;
        int step = 0;
        visited[1] = true;

        int answer = 0;

        while (true) {
            cur = A[cur];
            step++;

            if (cur == N) {
                answer = step;
                break;
            }

            if (visited[cur]) {
                answer = 0;
                break;
            }

            visited[cur] = true;
        }

        cout << answer << '\n';
    }

    return 0;
}