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
                dp[i + 1][rem] = true;
                pick[i + 1][rem] = false;
                int nxt = (rem + val) % (M - 1);
                dp[i + 1][nxt] = true;
                pick[i + 1][nxt] = true;
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
                dp[i + 1][rem] = true;
                pick[i + 1][rem] = false;
                int nxt = (rem + val) % (N - 1);
                dp[i + 1][nxt] = true;
                pick[i + 1][nxt] = true;
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

struct CompAssignment {
    long long delta;
    vector<int> choices;
};

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
    
    int V = (N - 1) + (M - 1);
    vector<vector<long long>> cand(V);
    
    for (int i = 1; i < N; ++i) {
        cand[i - 1].push_back(C[0][0] - C[i][0]);
        if (C[i][0] != 0) cand[i - 1].push_back(C[0][0] + C[i][0]);
    }
    for (int j = 1; j < M; ++j) {
        cand[N - 1 + j - 1].push_back(-C[0][j]);
        if (C[0][j] != 0) cand[N - 1 + j - 1].push_back(C[0][j]);
    }
    
    vector<vector<bool>> valid(V);
    for (int i = 0; i < V; ++i) {
        for (int opt = 0; opt < cand[i].size(); ++opt) valid[i].push_back(true);
    }
    
    // 1. 호 일관성(Arc Consistency) 적용: 절대 불가능한 선택지 제거
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 1; i < N; ++i) {
            int u = i - 1;
            for (int a = 0; a < cand[u].size(); ++a) {
                if (!valid[u][a]) continue;
                for (int j = 1; j < M; ++j) {
                    int v = N - 1 + j - 1;
                    bool match = false;
                    for (int b = 0; b < cand[v].size(); ++b) {
                        if (valid[v][b] && abs(cand[v][b] - cand[u][a]) == C[i][j]) { match = true; break; }
                    }
                    if (!match) { valid[u][a] = false; changed = true; break; }
                }
            }
        }
        for (int j = 1; j < M; ++j) {
            int v = N - 1 + j - 1;
            for (int b = 0; b < cand[v].size(); ++b) {
                if (!valid[v][b]) continue;
                for (int i = 1; i < N; ++i) {
                    int u = i - 1;
                    bool match = false;
                    for (int a = 0; a < cand[u].size(); ++a) {
                        if (valid[u][a] && abs(cand[v][b] - cand[u][a]) == C[i][j]) { match = true; break; }
                    }
                    if (!match) { valid[v][b] = false; changed = true; break; }
                }
            }
        }
    }
    
    // 2. 종속성 그래프 구축 (선택을 강제하는 관계)
    vector<vector<int>> adj(V);
    for (int i = 1; i < N; ++i) {
        int u = i - 1;
        for (int j = 1; j < M; ++j) {
            int v = N - 1 + j - 1;
            bool dependent = false;
            
            for (int a = 0; a < cand[u].size(); ++a) {
                if (!valid[u][a]) continue;
                int count = 0;
                for (int b = 0; b < cand[v].size(); ++b) {
                    if (valid[v][b] && abs(cand[v][b] - cand[u][a]) == C[i][j]) count++;
                }
                if (count == 1) dependent = true;
            }
            
            for (int b = 0; b < cand[v].size(); ++b) {
                if (!valid[v][b]) continue;
                int count = 0;
                for (int a = 0; a < cand[u].size(); ++a) {
                    if (valid[u][a] && abs(cand[v][b] - cand[u][a]) == C[i][j]) count++;
                }
                if (count == 1) dependent = true;
            }
            
            if (dependent) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
    }
    
    // 3. 진정한 의미의 독립 컴포넌트 추출
    vector<vector<int>> comps;
    vector<bool> vis(V, false);
    for (int i = 0; i < V; ++i) {
        if (!vis[i]) {
            vector<int> comp;
            queue<int> q;
            q.push(i); vis[i] = true;
            while (!q.empty()) {
                int curr = q.front(); q.pop();
                comp.push_back(curr);
                for (int nxt : adj[curr]) {
                    if (!vis[nxt]) { vis[nxt] = true; q.push(nxt); }
                }
            }
            comps.push_back(comp);
        }
    }
    
    // 4. 각 컴포넌트별 백트래킹으로 모든 유효한 선택지 수집
    vector<vector<CompAssignment>> comp_valid_assigns(comps.size());
    for (int c_idx = 0; c_idx < comps.size(); ++c_idx) {
        vector<int> current_choices;
        
        auto backtrack = [&](auto& self, int var_idx) -> void {
            if (var_idx == comps[c_idx].size()) {
                long long sR = 0, sC = 0;
                for (int i = 0; i < comps[c_idx].size(); ++i) {
                    int u = comps[c_idx][i];
                    if (u < N - 1) sR += cand[u][current_choices[i]];
                    else sC += cand[u][current_choices[i]];
                }
                comp_valid_assigns[c_idx].push_back({sC - sR, current_choices});
                return;
            }
            int u = comps[c_idx][var_idx];
            for (int opt = 0; opt < cand[u].size(); ++opt) {
                if (!valid[u][opt]) continue;
                bool ok = true;
                for (int i = 0; i < var_idx; ++i) {
                    int v = comps[c_idx][i];
                    if (u < N - 1 && v >= N - 1) {
                        if (abs(cand[v][current_choices[i]] - cand[u][opt]) != C[u + 1][v - (N - 1) + 1]) { ok = false; break; }
                    } else if (u >= N - 1 && v < N - 1) {
                        if (abs(cand[u][opt] - cand[v][current_choices[i]]) != C[v + 1][u - (N - 1) + 1]) { ok = false; break; }
                    }
                }
                if (ok) {
                    current_choices.push_back(opt);
                    self(self, var_idx + 1);
                    current_choices.pop_back();
                }
            }
        };
        backtrack(backtrack, 0);
    }
    
    // 5. 배낭 문제(DP)를 활용한 최적 조합 구성
    long long current_diff = C[0][0]; 
    vector<int> final_choices(V);
    
    if (N != M) {
        int mod = abs(N - M);
        vector<vector<bool>> dp(comps.size() + 1, vector<bool>(mod, false));
        vector<vector<int>> parent_choice(comps.size() + 1, vector<int>(mod, -1));
        
        int start_val = ((current_diff % mod) + mod) % mod;
        dp[0][start_val] = true;
        
        for (size_t i = 0; i < comps.size(); ++i) {
            for (int r = 0; r < mod; ++r) {
                if (!dp[i][r]) continue;
                for (size_t k = 0; k < comp_valid_assigns[i].size(); ++k) {
                    long long d = comp_valid_assigns[i][k].delta;
                    int nxt_r = ((r + d) % mod + mod) % mod;
                    dp[i + 1][nxt_r] = true;
                    parent_choice[i + 1][nxt_r] = k;
                }
            }
        }
        
        int curr_r = 0;
        for (int i = comps.size(); i > 0; --i) {
            int k = parent_choice[i][curr_r];
            auto& assign = comp_valid_assigns[i - 1][k];
            for (size_t j = 0; j < comps[i - 1].size(); ++j) {
                final_choices[comps[i - 1][j]] = assign.choices[j];
            }
            curr_r = ((curr_r - assign.delta) % mod + mod) % mod;
        }
    } else {
        vector<vector<pair<long long, int>>> dp(comps.size() + 1);
        dp[0].push_back({current_diff, -1});
        
        for (size_t i = 0; i < comps.size(); ++i) {
            unordered_map<long long, int> next_dp;
            for (auto& state : dp[i]) {
                for (size_t k = 0; k < comp_valid_assigns[i].size(); ++k) {
                    long long nxt_val = state.first + comp_valid_assigns[i][k].delta;
                    next_dp[nxt_val] = k;
                }
            }
            for (auto& kv : next_dp) dp[i + 1].push_back({kv.first, kv.second});
        }
        
        long long curr_val = 0;
        for (int i = comps.size(); i > 0; --i) {
            int k = -1;
            for (auto& state : dp[i]) {
                if (state.first == curr_val) { k = state.second; break; }
            }
            auto& assign = comp_valid_assigns[i - 1][k];
            for (size_t j = 0; j < comps[i - 1].size(); ++j) {
                final_choices[comps[i - 1][j]] = assign.choices[j];
            }
            curr_val -= assign.delta;
        }
    }
    
    // 6. 결과 복원 및 출력
    vector<long long> R(N), C_arr(M);
    R[0] = 0; C_arr[0] = C[0][0];
    for (int i = 1; i < N; ++i) R[i] = cand[i - 1][final_choices[i - 1]];
    for (int j = 1; j < M; ++j) C_arr[j] = cand[N - 1 + j - 1][final_choices[N - 1 + j - 1]];
    
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