#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MAXN = 100005;
const int MAXM = 1000005;

int head[MAXN];
int to[MAXM];
int nxt[MAXM];
int in_degree[MAXN];
int topo[MAXN];
int L[MAXN];
int R[MAXN];
int freqL[MAXN];
int freqR[MAXN];
int CountA[MAXN];
int CountD[MAXN];
uint64_t dp[MAXN][16];

int edge_cnt = 0;

void add_edge(int u, int v) {
    to[edge_cnt] = v;
    nxt[edge_cnt] = head[u];
    head[u] = edge_cnt++;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M_edges;
    if (!(cin >> N >> M_edges)) return 0;

    for (int i = 1; i <= N; ++i) {
        head[i] = -1;
    }

    for (int i = 0; i < M_edges; ++i) {
        int u, v;
        cin >> u >> v;
        add_edge(u, v);
        in_degree[v]++;
    }

    int q_head = 0, q_tail = 0;
    for (int i = 1; i <= N; ++i) {
        if (in_degree[i] == 0) {
            topo[q_tail++] = i;
        }
    }

    while (q_head < q_tail) {
        int u = topo[q_head++];
        for (int e = head[u]; e != -1; e = nxt[e]) {
            int v = to[e];
            if (--in_degree[v] == 0) {
                topo[q_tail++] = v;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        int u = topo[i];
        for (int e = head[u]; e != -1; e = nxt[e]) {
            int v = to[e];
            if (L[v] < L[u] + 1) L[v] = L[u] + 1;
        }
    }

    for (int i = N - 1; i >= 0; --i) {
        int u = topo[i];
        for (int e = head[u]; e != -1; e = nxt[e]) {
            int v = to[e];
            if (R[u] < R[v] + 1) R[u] = R[v] + 1;
        }
    }

    int max_len = 0;
    for (int i = 1; i <= N; ++i) {
        if (L[i] + R[i] > max_len) {
            max_len = L[i] + R[i];
        }
        freqL[L[i]]++;
        freqR[R[i]]++;
    }

    vector<int> cands;
    for (int i = 1; i <= N; ++i) {
        if (L[i] + R[i] == max_len && freqL[L[i]] == 1 && freqR[R[i]] == 1) {
            cands.push_back(i);
        }
    }

    if (cands.empty()) {
        cout << 0 << "\n";
        return 0;
    }

    int num_cands = cands.size();
    int BATCH_BITS = 1024;
    int BATCH_WORDS = 16;
    int num_batches = (num_cands + BATCH_BITS - 1) / BATCH_BITS;

    for (int batch = 0; batch < num_batches; ++batch) {
        int start_idx = batch * BATCH_BITS;
        int end_idx = min(num_cands, start_idx + BATCH_BITS);

        for (int i = 1; i <= N; ++i) {
            for (int w = 0; w < BATCH_WORDS; ++w) {
                dp[i][w] = 0;
            }
        }
        for (int i = start_idx; i < end_idx; ++i) {
            int bit_idx = i - start_idx;
            dp[cands[i]][bit_idx >> 6] |= (1ULL << (bit_idx & 63));
        }

        for (int i = N - 1; i >= 0; --i) {
            int u = topo[i];
            for (int e = head[u]; e != -1; e = nxt[e]) {
                int v = to[e];
                for (int w = 0; w < BATCH_WORDS; ++w) {
                    dp[u][w] |= dp[v][w];
                }
            }
        }

        for (int v = 1; v <= N; ++v) {
            for (int w = 0; w < BATCH_WORDS; ++w) {
                uint64_t mask = dp[v][w];
                while (mask) {
                    int bit = __builtin_ctzll(mask);
                    CountA[start_idx + w * 64 + bit]++;
                    mask &= mask - 1;
                }
            }
        }

        for (int i = 1; i <= N; ++i) {
            for (int w = 0; w < BATCH_WORDS; ++w) {
                dp[i][w] = 0;
            }
        }
        for (int i = start_idx; i < end_idx; ++i) {
            int bit_idx = i - start_idx;
            dp[cands[i]][bit_idx >> 6] |= (1ULL << (bit_idx & 63));
        }

        for (int i = 0; i < N; ++i) {
            int u = topo[i];
            for (int e = head[u]; e != -1; e = nxt[e]) {
                int v = to[e];
                for (int w = 0; w < BATCH_WORDS; ++w) {
                    dp[v][w] |= dp[u][w];
                }
            }
        }

        for (int v = 1; v <= N; ++v) {
            for (int w = 0; w < BATCH_WORDS; ++w) {
                uint64_t mask = dp[v][w];
                while (mask) {
                    int bit = __builtin_ctzll(mask);
                    CountD[start_idx + w * 64 + bit]++;
                    mask &= mask - 1;
                }
            }
        }
    }

    vector<int> ans;
    for (int i = 0; i < num_cands; ++i) {
        if (CountA[i] + CountD[i] - 1 == N) {
            ans.push_back(cands[i]);
        }
    }

    sort(ans.begin(), ans.end());
    
    cout << ans.size() << "\n";
    if (ans.empty()) {
        cout << 0 << "\n";
    } else {
        for (int i = 0; i < ans.size(); ++i) {
            cout << ans[i] << (i + 1 == ans.size() ? "" : " ");
        }
        cout << "\n";
    }

    return 0;
}