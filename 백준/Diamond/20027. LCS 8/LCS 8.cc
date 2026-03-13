#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

const int MOD = 1e9 + 7;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string S;
    if (!(cin >> S)) return 0;
    int K;
    cin >> K;

    int N = S.length();

    map<vector<int>, int> state_to_id;
    vector<vector<int>> id_to_state;
    vector<vector<int>> tr;

    auto get_id = [&](vector<int> V) {
        if (state_to_id.count(V)) return state_to_id[V];
        int id = id_to_state.size();
        state_to_id[V] = id;
        id_to_state.push_back(V);
        tr.push_back(vector<int>(1 << (2 * K + 1), -1));
        return id;
    };

    vector<int> init_V(2 * K + 1, -100);
    for (int d = -K; d <= K; ++d) {
        int L = (d >= 0) ? 0 : -100;
        int min_req = max(0, d) - K;
        if (L >= min_req) {
            init_V[d + K] = L;
        }
    }
    get_id(init_V);

    int head = 0;
    while (head < id_to_state.size()) {
        vector<int> V = id_to_state[head];
        for (int mask = 0; mask < (1 << (2 * K + 1)); ++mask) {
            vector<int> next_V(2 * K + 1, -100);
            bool valid = false;
            for (int d = -K; d <= K; ++d) {
                int v1 = -100, v2 = -100, v3 = -100;
                if (d + 1 <= K) {
                    v1 = (V[d + 1 + K] == -100) ? -100 : V[d + 1 + K] - 1;
                }
                if (d - 1 >= -K) {
                    v2 = next_V[d - 1 + K];
                }
                if ((mask >> (d + K)) & 1) {
                    v3 = V[d + K];
                }
                int mx = max({v1, v2, v3});
                int min_req = max(0, d) - K;
                if (mx >= min_req) {
                    next_V[d + K] = min(mx, 0);
                    valid = true;
                }
            }
            if (valid) {
                tr[head][mask] = get_id(next_V);
            }
        }
        head++;
    }

    int num_states = id_to_state.size();
    vector<int> dp(num_states, 0);
    dp[0] = 1;
    vector<int> next_dp(num_states, 0);

    for (int i = 0; i < N; ++i) {
        int masks[26] = {0};
        for (int c = 0; c < 26; ++c) {
            int mask = 0;
            for (int d = -K; d <= K; ++d) {
                int j = i + d;
                if (j >= 0 && j < N && S[j] == c + 'A') {
                    mask |= (1 << (d + K));
                }
            }
            masks[c] = mask;
        }

        fill(next_dp.begin(), next_dp.end(), 0);
        for (int u = 0; u < num_states; ++u) {
            if (!dp[u]) continue;
            for (int c = 0; c < 26; ++c) {
                int v = tr[u][masks[c]];
                if (v != -1) {
                    next_dp[v] += dp[u];
                    if (next_dp[v] >= MOD) next_dp[v] -= MOD;
                }
            }
        }
        dp = next_dp;
    }

    int ans = 0;
    for (int u = 0; u < num_states; ++u) {
        if (dp[u] > 0 && id_to_state[u][K] != -100) {
            ans = (ans + dp[u]) % MOD;
        }
    }

    cout << ans << "\n";

    return 0;
}