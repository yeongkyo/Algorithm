#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
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
            if (dp[M][r]) {
                target_rem = r;
                break;
            }
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
    
    for (int j = 0; j < M; ++j) {
        cout << s[j] * C[0][j] + X << (j == M - 1 ? "" : " ");
    }
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
            if (dp[N][r]) {
                target_rem = r;
                break;
            }
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
    
    for (int i = 0; i < N; ++i) {
        cout << s[i] * C[i][0] + X << "\n";
    }
}

struct Implication {
    int v, b;
};

vector<Implication> implies_graph[2005][2];
bool cand_valid[2005][2];

bool dfs_assign(int x, int c, vector<int>& local_assign) {
    if (local_assign[x] != -1) return local_assign[x] == c;
    if (!cand_valid[x][c]) return false;
    local_assign[x] = c;
    for (auto& imp : implies_graph[x][c]) {
        if (!dfs_assign(imp.v, imp.b, local_assign)) return false;
    }
    return true;
}

struct Component {
    long long sumR0, sumC0, sumR1, sumC1;
    vector<pair<int, int>> assigns0, assigns1;
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
    
    int num_vars = (N - 1) + (M - 1);
    vector<vector<long long>> cand(num_vars);
    
    for (int i = 1; i < N; ++i) {
        cand[i - 1].push_back(C[0][0] - C[i][0]);
        if (C[i][0] != 0) cand[i - 1].push_back(C[0][0] + C[i][0]);
    }
    for (int j = 1; j < M; ++j) {
        cand[N - 1 + j - 1].push_back(-C[0][j]);
        if (C[0][j] != 0) cand[N - 1 + j - 1].push_back(C[0][j]);
    }
    
    for (int i = 0; i < num_vars; ++i) {
        cand_valid[i][0] = true;
        cand_valid[i][1] = (cand[i].size() > 1);
    }
    
    for (int i = 1; i < N; ++i) {
        int u = i - 1;
        for (int j = 1; j < M; ++j) {
            int v = N - 1 + j - 1;
            
            for (int a = 0; a < cand[u].size(); ++a) {
                int valid_b = -1, count = 0;
                for (int b = 0; b < cand[v].size(); ++b) {
                    if (abs(cand[v][b] - cand[u][a]) == C[i][j]) { valid_b = b; count++; }
                }
                if (count == 0) cand_valid[u][a] = false;
                else if (count == 1) implies_graph[u][a].push_back({v, valid_b});
            }
            
            for (int b = 0; b < cand[v].size(); ++b) {
                int valid_a = -1, count = 0;
                for (int a = 0; a < cand[u].size(); ++a) {
                    if (abs(cand[v][b] - cand[u][a]) == C[i][j]) { valid_a = a; count++; }
                }
                if (count == 0) cand_valid[v][b] = false;
                else if (count == 1) implies_graph[v][b].push_back({u, valid_a});
            }
        }
    }
    
    vector<int> global_assign(num_vars, -1);
    
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < num_vars; ++i) {
            if (global_assign[i] != -1) continue;
            vector<int> a0 = global_assign; bool ok0 = dfs_assign(i, 0, a0);
            vector<int> a1 = global_assign; bool ok1 = false;
            if (cand[i].size() > 1) ok1 = dfs_assign(i, 1, a1);
            
            if (ok0 && !ok1) { global_assign = a0; changed = true; }
            else if (!ok0 && ok1) { global_assign = a1; changed = true; }
        }
    }
    
    vector<int> forced_assign = global_assign;
    long long base_sumR = 0, base_sumC = C[0][0];
    vector<int> final_assign(num_vars, -1);
    
    for (int i = 0; i < num_vars; ++i) {
        if (forced_assign[i] != -1) {
            final_assign[i] = forced_assign[i];
            if (i < N - 1) base_sumR += cand[i][forced_assign[i]];
            else base_sumC += cand[i][forced_assign[i]];
        }
    }
    
    vector<Component> components;
    for (int i = 0; i < num_vars; ++i) {
        if (global_assign[i] != -1) continue;
        vector<int> a0 = global_assign; bool ok0 = dfs_assign(i, 0, a0);
        vector<int> a1 = global_assign; bool ok1 = dfs_assign(i, 1, a1);
        
        Component comp;
        comp.sumR0 = comp.sumC0 = comp.sumR1 = comp.sumC1 = 0;
        
        for (int j = 0; j < num_vars; ++j) {
            if (global_assign[j] == -1) {
                if (a0[j] != -1) {
                    comp.assigns0.push_back({j, a0[j]});
                    if (j < N - 1) comp.sumR0 += cand[j][a0[j]];
                    else comp.sumC0 += cand[j][a0[j]];
                }
                if (a1[j] != -1) {
                    comp.assigns1.push_back({j, a1[j]});
                    if (j < N - 1) comp.sumR1 += cand[j][a1[j]];
                    else comp.sumC1 += cand[j][a1[j]];
                }
            }
        }
        components.push_back(comp);
        for (auto& p : comp.assigns0) global_assign[p.first] = p.second;
    }
    
    long long current_diff = base_sumC - base_sumR;
    vector<long long> deltas;
    for (auto& comp : components) {
        long long diff0 = comp.sumC0 - comp.sumR0;
        long long diff1 = comp.sumC1 - comp.sumR1;
        current_diff += diff0;
        deltas.push_back(diff1 - diff0);
    }
    
    if (N != M) {
        int mod = abs(N - M);
        vector<vector<bool>> dp(components.size() + 1, vector<bool>(mod, false));
        vector<vector<bool>> pick(components.size() + 1, vector<bool>(mod, false));
        
        int start_val = ((current_diff % mod) + mod) % mod;
        dp[0][start_val] = true;
        
        for (size_t i = 0; i < deltas.size(); ++i) {
            long long d = ((deltas[i] % mod) + mod) % mod;
            for (int r = 0; r < mod; ++r) {
                if (dp[i][r]) {
                    dp[i + 1][r] = true;
                    pick[i + 1][r] = false;
                    dp[i + 1][(r + d) % mod] = true;
                    pick[i + 1][(r + d) % mod] = true;
                }
            }
        }
        
        int curr_r = 0;
        for (int i = components.size(); i > 0; --i) {
            if (pick[i][curr_r]) {
                for (auto& p : components[i - 1].assigns1) final_assign[p.first] = p.second;
                curr_r = (curr_r - (((deltas[i - 1] % mod) + mod) % mod) + mod) % mod;
            } else {
                for (auto& p : components[i - 1].assigns0) final_assign[p.first] = p.second;
            }
        }
    } else {
        vector<unordered_map<long long, int>> dp_exact(components.size() + 1);
        dp_exact[0][current_diff] = 0;
        
        for (size_t i = 0; i < deltas.size(); ++i) {
            long long d = deltas[i];
            for (auto& kv : dp_exact[i]) dp_exact[i + 1][kv.first] = 0;
            for (auto& kv : dp_exact[i]) dp_exact[i + 1][kv.first + d] = 1;
        }
        
        long long curr_val = 0;
        for (int i = components.size(); i > 0; --i) {
            if (dp_exact[i][curr_val] == 1) {
                for (auto& p : components[i - 1].assigns1) final_assign[p.first] = p.second;
                curr_val -= deltas[i - 1];
            } else {
                for (auto& p : components[i - 1].assigns0) final_assign[p.first] = p.second;
            }
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
    for(int i = 0; i < N - 1; ++i) A[i][M - 1] = R[i];
    for(int j = 0; j < M - 1; ++j) A[N - 1][j] = C_arr[j];
    A[N - 1][M - 1] = C_arr[M - 1] + R[N - 1] - sumC;
    
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j) {
            cout << A[i][j] << (j == M - 1 ? "" : " ");
        }
        cout << "\n";
    }
    
    return 0;
}