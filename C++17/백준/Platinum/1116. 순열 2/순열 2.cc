#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> P(N);
    vector<int> prev(N);
    for (int i = 0; i < N; ++i) {
        cin >> P[i];
        prev[P[i]] = i;
    }

    vector<bool> visited(N, false);
    vector<int> C0;
    vector<int> m; 

    for (int i = 0; i < N; ++i) {
        if (!visited[i]) {
            vector<int> cycle;
            int curr = i;
            while (!visited[curr]) {
                visited[curr] = true;
                cycle.push_back(curr);
                curr = P[curr];
            }
            if (find(cycle.begin(), cycle.end(), 0) != cycle.end()) {
                C0 = cycle;
            } else {
                int min_val = cycle[0];
                for (int val : cycle) min_val = min(min_val, val);
                m.push_back(min_val);
            }
        }
    }

    if (m.empty()) {
        for (int i = 0; i < N; ++i) {
            cout << P[i] << (i == N - 1 ? "" : " ");
        }
        cout << "\n";
        return 0;
    }

    sort(m.begin(), m.end());

    vector<int> best_B(N, 1e9);
    vector<int> best_Q;

    for (int x : C0) {
        vector<int> B;
        int curr = 0;
        while (true) {
            B.push_back(curr);
            if (curr == x) break;
            curr = P[curr];
        }

        for (int i = 0; i < m.size(); ++i) {
            curr = m[i];
            while (true) {
                B.push_back(curr);
                if (curr == prev[m[i]]) break;
                curr = P[curr];
            }
        }

        curr = P[x];
        while (B.size() < N) {
            B.push_back(curr);
            curr = P[curr];
        }

        if (B < best_B) {
            best_B = B;
            vector<int> Q = P;
            Q[x] = m[0];
            for (int i = 0; i < (int)m.size() - 1; ++i) {
                Q[prev[m[i]]] = m[i + 1];
            }
            Q[prev[m.back()]] = P[x];
            best_Q = Q;
        }
    }

    for (int i = 0; i < N; ++i) {
        cout << best_Q[i] << (i == N - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}