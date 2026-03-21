#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const long long INF = 2e18;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long M;
    if (!(cin >> M)) return 0;

    int N;
    cin >> N;

    vector<long long> A(N);
    long long K = 0;
    for (int i = 0; i < N; ++i) {
        cin >> A[i];
        if (A[i] > K) {
            K = A[i];
        }
    }

    sort(A.begin(), A.end());
    A.erase(unique(A.begin(), A.end()), A.end());

    vector<long long> dist(K, INF);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;

    dist[0] = 0;
    pq.push({0, 0});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        for (long long c : A) {
            int nxt = (u + c) % K;
            long long cost = K - c;
            
            if (dist[u] + cost < dist[nxt]) {
                dist[nxt] = dist[u] + cost;
                pq.push({dist[nxt], nxt});
            }
        }
    }

    long long r = M % K;
    long long ans = (M + dist[r]) / K;

    cout << ans << "\n";

    return 0;
}