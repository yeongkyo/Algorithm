#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <unordered_map>

using namespace std;

int N, M;
long long C[1005][1005];

void solve_1_M() {
    long long S = 0;
    for (int j = 0; j < M; ++j) S += C[0][j];
    
    vector<vector<bool>> dp(M + 1, vector<bool>(M - 1, false));
    vector<vector<bool>> pick(M + 1, vector<bool>(M - 1, false));
    dp[0][0] = true;
    
    for (int i = 0; i < M; ++i) {
        int val = C[0][i] % (M - 1);
        for (int rem = 0; rem < M - 1; ++rem) {
            if (dp[i][rem]) {
                dp[i + 1][rem] = true; pick[i + 1][rem] = false;
                dp[i + 1][(rem + val) % (M - 1)] = true; pick[i + 1][(rem + val) % (M - 1)] = true;
            }
        }
    }
    
    int target_rem = -1;
    for (int r = 0; r < M - 1; ++r) {
        if ((2 * r) % (M - 1) == S % (M - 1)) {
            if (dp[M][r]) { target_rem = r; break; }
        }
    }
    
    vector<int> s(M, 1);
    int curr = target_rem;
    for (int i = M; i > 0; --i) {
        if (pick[i][curr]) {
            s[i - 1] = -1;
            curr = (curr - (C[0][i - 1] % (M - 1)) + (M - 1)) % (M - 1);
        }
    }
    long long sum_sC = 0;
    for (int j = 0; j < M; ++j) sum_sC += s[j] * C[0][j];
    long long X = sum_sC / (1 - M);
    for (int j = 0; j < M; ++j) cout << s[j] * C[0][j] + X << (j == M - 1 ? "" : " ");
    cout << "\n";
}

void solve_N_1() {
    long long S = 0;
    for (int i = 0; i < N; ++i) S += C[i][0];
    
    vector<vector<bool>> dp(N + 1, vector<bool>(N - 1, false));
    vector<vector<bool>> pick(N + 1, vector<bool>(N - 1, false));
    dp[0][0] = true;
    
    for (int i = 0; i < N; ++i) {
        int val = C[i][0] % (N - 1);
        for (int rem = 0; rem < N - 1; ++rem) {
            if (dp[i][rem]) {
                dp[i + 1][rem] = true; pick[i + 1][rem] = false;
                dp[i + 1][(rem + val) % (N - 1)] = true; pick[i + 1][(rem + val) % (N - 1)] = true;
            }
        }
    }
    
    int target_rem = -1;
    for (int r = 0; r < N - 1; ++r) {
        if ((2 * r) % (N - 1) == S % (N - 1)) {
            if (dp[N][r]) { target_rem = r; break; }
        }
    }
    
    vector<int> s(N, 1);
    int curr = target_rem;
    for (int i = N; i > 0; --i) {
        if (pick[i][curr]) {
            s[i - 1] = -1;
            curr = (curr - (C[i - 1][0] % (N - 1)) + (N - 1)) % (N - 1);
        }
    }
    long long sum_sC = 0;
    for (int i = 0; i < N; ++i) sum_sC += s[i] * C[i][0];
    long long X = sum_sC / (1 - N);
    for (int i = 0; i < N; ++i) cout << s[i] * C[i][0] + X << "\n";
}

struct MacroState {
    long long base_delta;
    vector<long long> free_deltas;
    vector<pair<int, int>> base_assign;
    vector<pair<int, int>> free_vars;
};

struct Transition {
    int prev_r;
    int ms_idx;
    vector<bool> free_choices;
};

struct TransitionExact {
    long long prev_r;
    int ms_idx;
    vector<bool> free_choices;
};

int V;
vector<vector<long long>> cand;
vector<vector<int>> incompatible[2000];
int assign_arr[2000];
int invalid_count[2000][2];
vector<MacroState> current_comp_states;

bool assign_val(int u, int a, vector<int>& trail) {
    assign_arr[u] = a;
    trail.push_back(u);
    for (int enc : incompatible[u][a]) {
        int v = enc / 2, b = enc % 2;
        invalid_count[v][b]++;
        if (assign_arr[v] == -1 && invalid_count[v][0] > 0 && invalid_count[v][1] > 0) return false;
    }
    for (int enc : incompatible[u][a]) {
        int v = enc / 2;
        if (assign_arr[v] == -1) {
            if (invalid_count[v][0] > 0) {
                if (!assign_val(v, 1, trail)) return false;
            } else if (invalid_count[v][1] > 0) {
                if (!assign_val(v, 0, trail)) return false;
            }
        }
    }
    return true;
}

void backtrack_comp(const vector<int>& unassigned) {
    int branch_var = -1;
    for (int u : unassigned) {
        if (assign_arr[u] != -1) continue;
        bool has_neighbor = false;
        for (int a = 0; a < 2; ++a) {
            if (invalid_count[u][a] > 0) continue;
            for (int enc : incompatible[u][a]) {
                if (assign_arr[enc / 2] == -1) { has_neighbor = true; break; }
            }
            if (has_neighbor) break;
        }
        if (has_neighbor) { branch_var = u; break; }
    }
    
    if (branch_var == -1) {
        MacroState state; state.base_delta = 0;
        for (int u : unassigned) {
            if (assign_arr[u] != -1) {
                long long val = cand[u][assign_arr[u]];
                state.base_delta += (u < N - 1) ? -val : val;
                state.base_assign.push_back({u, assign_arr[u]});
            } else {
                long long d0 = (u < N - 1) ? -cand[u][0] : cand[u][0];
                long long d1 = (u < N - 1) ? -cand[u][1] : cand[u][1];
                state.base_delta += d0;
                state.free_deltas.push_back(d1 - d0);
                state.base_assign.push_back({u, 0});
                state.free_vars.push_back({u, 1});
            }
        }
        current_comp_states.push_back(state);
        return;
    }
    
    for (int a = 0; a < 2; ++a) {
        if (invalid_count[branch_var][a] == 0) {
            vector<int> trail;
            if (assign_val(branch_var, a, trail)) backtrack_comp(unassigned);
            for (int v : trail) {
                int val = assign_arr[v];
                assign_arr[v] = -1;
                for (int enc : incompatible[v][val]) invalid_count[enc / 2][enc % 2]--;
            }
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N >> M)) return 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) cin >> C[i][j];
    }
    
    if (N == 1 && M == 1) { cout << 0 << "\n"; return 0; }
    if (N == 1) { solve_1_M(); return 0; }
    if (M == 1) { solve_N_1(); return 0; }
    
    V = (N - 1) + (M - 1);
    cand.resize(V);
    for (int i = 1; i < N; ++i) {
        cand[i - 1].push_back(C[0][0] - C[i][0]);
        cand[i - 1].push_back(C[0][0] + C[i][0]);
    }
    for (int j = 1; j < M; ++j) {
        cand[N - 1 + j - 1].push_back(-C[0][j]);
        cand[N - 1 + j - 1].push_back(C[0][j]);
    }
    
    for (int i = 0; i < V; ++i) incompatible[i].resize(2);
    
    for (int i = 1; i < N; ++i) {
        int u = i - 1;
        for (int j = 1; j < M; ++j) {
            int v = N - 1 + j - 1;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    if (abs(cand[v][b] - cand[u][a]) != C[i][j]) {
                        incompatible[u][a].push_back(v * 2 + b);
                        incompatible[v][b].push_back(u * 2 + a);
                    }
                }
            }
        }
    }
    
    fill(assign_arr, assign_arr + V, -1);
    bool initial_valid[2000][2];
    for (int i = 0; i < V; ++i) { initial_valid[i][0] = true; initial_valid[i][1] = true; }
    bool changed = true;
    while (changed) {
        changed = false;
        for (int u = 0; u < V; ++u) {
            for (int a = 0; a < 2; ++a) {
                if (!initial_valid[u][a]) continue;
                for (int v = 0; v < V; ++v) {
                    if (u == v) continue;
                    if ((u < N - 1 && v >= N - 1) || (u >= N - 1 && v < N - 1)) {
                        bool can_match = false;
                        for (int b = 0; b < 2; ++b) {
                            if (initial_valid[v][b]) {
                                bool incompat = false;
                                for (int enc : incompatible[u][a]) { if (enc == v * 2 + b) incompat = true; }
                                if (!incompat) can_match = true;
                            }
                        }
                        if (!can_match) { initial_valid[u][a] = false; changed = true; break; }
                    }
                }
            }
        }
    }
    
    for (int u = 0; u < V; ++u) {
        if (!initial_valid[u][0]) { vector<int> t; assign_val(u, 1, t); }
        else if (!initial_valid[u][1]) { vector<int> t; assign_val(u, 0, t); }
    }
    
    vector<vector<MacroState>> all_comp_states;
    vector<bool> vis(V, false);
    for (int i = 0; i < V; ++i) {
        if (vis[i] || assign_arr[i] != -1) continue;
        vector<int> comp; queue<int> q;
        q.push(i); vis[i] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int a = 0; a < 2; ++a) {
                if (invalid_count[u][a] > 0) continue;
                for (int enc : incompatible[u][a]) {
                    int v = enc / 2;
                    if (!vis[v] && assign_arr[v] == -1) { vis[v] = true; q.push(v); }
                }
            }
        }
        current_comp_states.clear();
        backtrack_comp(comp);
        all_comp_states.push_back(current_comp_states);
    }
    
    long long global_base_sumR = 0, global_base_sumC = C[0][0];
    vector<int> final_assign(V, -1);
    for (int i = 0; i < V; ++i) {
        if (assign_arr[i] != -1) {
            final_assign[i] = assign_arr[i];
            if (i < N - 1) global_base_sumR += cand[i][assign_arr[i]];
            else global_base_sumC += cand[i][assign_arr[i]];
        }
    }
    
    long long start_delta = global_base_sumC - global_base_sumR;
    
    if (N != M) {
        int mod = abs(N - M);
        vector<vector<Transition>> dp_hist(all_comp_states.size() + 1, vector<Transition>(mod, {-1, -1, {}}));
        int start_mod = ((start_delta % mod) + mod) % mod;
        dp_hist[0][start_mod].prev_r = -2;
        
        for (size_t c = 0; c < all_comp_states.size(); ++c) {
            for (int ms_idx = 0; ms_idx < all_comp_states[c].size(); ++ms_idx) {
                auto& ms = all_comp_states[c][ms_idx];
                int K = ms.free_deltas.size();
                vector<vector<int>> track_choice(K + 1, vector<int>(mod, -1));
                vector<vector<int>> track_prev(K + 1, vector<int>(mod, -1));
                
                for (int r = 0; r < mod; ++r) {
                    if (dp_hist[c][r].prev_r != -1) {
                        int start_r = ((r + ms.base_delta) % mod + mod) % mod;
                        track_prev[0][start_r] = r;
                    }
                }
                for (int k = 0; k < K; ++k) {
                    long long d = ((ms.free_deltas[k] % mod) + mod) % mod;
                    for (int r = 0; r < mod; ++r) {
                        if (track_prev[k][r] != -1) {
                            track_prev[k+1][r] = r; track_choice[k+1][r] = 0;
                        }
                    }
                    for (int r = 0; r < mod; ++r) {
                        if (track_prev[k][r] != -1) {
                            int nxt_r = (r + d) % mod;
                            track_prev[k+1][nxt_r] = r; track_choice[k+1][nxt_r] = 1;
                        }
                    }
                }
                for (int r = 0; r < mod; ++r) {
                    if (track_prev[K][r] != -1 && dp_hist[c+1][r].prev_r == -1) {
                        Transition t; t.ms_idx = ms_idx; t.free_choices.resize(K);
                        int curr = r;
                        for (int k = K; k > 0; --k) {
                            t.free_choices[k-1] = track_choice[k][curr];
                            curr = track_prev[k][curr];
                        }
                        t.prev_r = track_prev[0][curr];
                        dp_hist[c+1][r] = t;
                    }
                }
            }
        }
        int curr = 0;
        for (int c = all_comp_states.size(); c > 0; --c) {
            auto& t = dp_hist[c][curr];
            auto& ms = all_comp_states[c - 1][t.ms_idx];
            for (auto& p : ms.base_assign) final_assign[p.first] = p.second;
            for (size_t k = 0; k < t.free_choices.size(); ++k) {
                if (t.free_choices[k]) final_assign[ms.free_vars[k].first] = ms.free_vars[k].second;
            }
            curr = t.prev_r;
        }
    } else {
        vector<unordered_map<long long, TransitionExact>> dp_hist(all_comp_states.size() + 1);
        dp_hist[0][start_delta] = {-2LL, -1, {}};
        
        for (size_t c = 0; c < all_comp_states.size(); ++c) {
            for (int ms_idx = 0; ms_idx < all_comp_states[c].size(); ++ms_idx) {
                auto& ms = all_comp_states[c][ms_idx];
                int K = ms.free_deltas.size();
                unordered_map<long long, long long> track_prev[K + 1];
                unordered_map<long long, int> track_choice[K + 1];
                
                for (auto& kv : dp_hist[c]) track_prev[0][kv.first + ms.base_delta] = kv.first;
                
                for (int k = 0; k < K; ++k) {
                    long long d = ms.free_deltas[k];
                    for (auto& kv : track_prev[k]) {
                        track_prev[k+1][kv.first] = kv.first; track_choice[k+1][kv.first] = 0;
                        track_prev[k+1][kv.first + d] = kv.first; track_choice[k+1][kv.first + d] = 1;
                    }
                }
                for (auto& kv : track_prev[K]) {
                    if (dp_hist[c+1].find(kv.first) == dp_hist[c+1].end()) {
                        TransitionExact t; t.ms_idx = ms_idx; t.free_choices.resize(K);
                        long long curr = kv.first;
                        for (int k = K; k > 0; --k) {
                            t.free_choices[k-1] = track_choice[k][curr];
                            curr = track_prev[k][curr];
                        }
                        t.prev_r = track_prev[0][curr];
                        dp_hist[c+1][kv.first] = t;
                    }
                }
            }
        }
        long long curr = 0;
        for (int c = all_comp_states.size(); c > 0; --c) {
            auto& t = dp_hist[c][curr];
            auto& ms = all_comp_states[c - 1][t.ms_idx];
            for (auto& p : ms.base_assign) final_assign[p.first] = p.second;
            for (size_t k = 0; k < t.free_choices.size(); ++k) {
                if (t.free_choices[k]) final_assign[ms.free_vars[k].first] = ms.free_vars[k].second;
            }
            curr = t.prev_r;
        }
    }
    
    vector<long long> R(N), C_arr(M);
    R[0] = 0; C_arr[0] = C[0][0];
    for (int i = 1; i < N; ++i) R[i] = cand[i - 1][final_assign[i - 1]];
    for (int j = 1; j < M; ++j) C_arr[j] = cand[N - 1 + j - 1][final_assign[N - 1 + j - 1]];
    
    long long sumR = 0, sumC = 0;
    for (long long x : R) sumR += x;
    for (long long x : C_arr) sumC += x;
    
    if (N != M) {
        long long Y = (sumC - sumR) / (N - M);
        for (int i = 0; i < N; ++i) R[i] += Y;
        for (int j = 0; j < M; ++j) C_arr[j] += Y;
        sumC += M * Y;
    }
    
    vector<vector<long long>> A(N, vector<long long>(M, 0));
    for (int i = 0; i < N - 1; ++i) A[i][M - 1] = R[i];
    for (int j = 0; j < M - 1; ++j) A[N - 1][j] = C_arr[j];
    A[N - 1][M - 1] = C_arr[M - 1] + R[N - 1] - sumC;
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) cout << A[i][j] << (j == M - 1 ? "" : " ");
        cout << "\n";
    }
    
    return 0;
}