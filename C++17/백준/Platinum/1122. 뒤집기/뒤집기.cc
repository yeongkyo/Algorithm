#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int A, B, K;
    if (!(cin >> A >> B >> K)) return 0;

    int N = A + B;
    vector<int> dist(N + 1, -1);
    set<int> unv[2];

    for (int i = 0; i <= N; ++i) {
        if (i != A) {
            unv[i % 2].insert(i);
        }
    }

    queue<int> q;
    dist[A] = 0;
    q.push(A);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == 0) break;

        int min_x = max(0, K - (N - u));
        int max_x = min(u, K);

        if (min_x > max_x) continue;

        int min_next = u + K - 2 * max_x;
        int max_next = u + K - 2 * min_x;
        int p = (u + K) % 2;

        auto it = unv[p].lower_bound(min_next);
        while (it != unv[p].end() && *it <= max_next) {
            int v = *it;
            dist[v] = dist[u] + 1;
            q.push(v);
            it = unv[p].erase(it);
        }
    }

    cout << dist[0] << "\n";

    return 0;
}