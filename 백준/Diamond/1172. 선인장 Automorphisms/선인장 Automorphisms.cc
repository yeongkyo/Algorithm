#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>
#include <map>
#include <cstring>

using namespace std;

long long MOD = 1000000003;
long long fact[405];

int N, M;
vector<vector<int>> adj;
int depth_arr[405], parent_arr[405];
bool visited[405];

vector<vector<int>> cycles;
vector<pair<int, int>> bridges;
set<pair<int, int>> cycle_edges;

int V_cnt, C_cnt, B_cnt;
vector<vector<int>> adj_T;

void dfs(int u, int p, int d) {
    visited[u] = true;
    depth_arr[u] = d;
    parent_arr[u] = p;
    for (int v : adj[u]) {
        if (v == p) continue;
        if (visited[v]) {
            if (depth_arr[v] < depth_arr[u]) {
                vector<int> cycle;
                int curr = u;
                while (curr != v) {
                    cycle.push_back(curr - 1);
                    curr = parent_arr[curr];
                }
                cycle.push_back(v - 1);
                cycles.push_back(cycle);
            }
        } else {
            dfs(v, u, d + 1);
        }
    }
}

int dist_T[1005], par_T[1005];
int bfs(int start) {
    memset(dist_T, -1, sizeof(dist_T));
    queue<int> q;
    q.push(start);
    dist_T[start] = 0;
    par_T[start] = -1;
    int farthest = start;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        farthest = u;
        for (int v : adj_T[u]) {
            if (dist_T[v] == -1) {
                dist_T[v] = dist_T[u] + 1;
                par_T[v] = u;
                q.push(v);
            }
        }
    }
    return farthest;
}

map<pair<int, vector<int>>, int> vec_to_id;
int get_id(int type, vector<int>& v) {
    auto key = make_pair(type, v);
    if (vec_to_id.find(key) == vec_to_id.end()) {
        vec_to_id[key] = vec_to_id.size() + 1;
    }
    return vec_to_id[key];
}

struct Result {
    int hash_id;
    long long ways;
};

Result HashWays(int u, int p) {
    if (u < N) {
        vector<int> c_hashes;
        long long ways = 1;
        for (int v : adj_T[u]) {
            if (v == p) continue;
            Result res = HashWays(v, u);
            c_hashes.push_back(res.hash_id);
            ways = (ways * res.ways) % MOD;
        }
        sort(c_hashes.begin(), c_hashes.end());
        int i = 0;
        while (i < c_hashes.size()) {
            int j = i;
            while (j < c_hashes.size() && c_hashes[j] == c_hashes[i]) j++;
            int count = j - i;
            ways = (ways * fact[count]) % MOD;
            i = j;
        }
        int my_hash = get_id(0, c_hashes);
        return {my_hash, ways};
    } else if (u >= N + C_cnt) {
        vector<int> c_hashes;
        long long ways = 1;
        for (int v : adj_T[u]) {
            if (v == p) continue;
            Result res = HashWays(v, u);
            c_hashes.push_back(res.hash_id);
            ways = (ways * res.ways) % MOD;
        }
        sort(c_hashes.begin(), c_hashes.end());
        int i = 0;
        while (i < c_hashes.size()) {
            int j = i;
            while (j < c_hashes.size() && c_hashes[j] == c_hashes[i]) j++;
            int count = j - i;
            ways = (ways * fact[count]) % MOD;
            i = j;
        }
        int my_hash = get_id(1, c_hashes);
        return {my_hash, ways};
    } else {
        int c_idx = u - N;
        vector<int>& cycle = cycles[c_idx];
        int k = cycle.size();
        
        if (p != -1) {
            int p_idx = -1;
            for (int i = 0; i < k; ++i) {
                if (cycle[i] == p) {
                    p_idx = i;
                    break;
                }
            }
            vector<int> seq_hashes;
            long long ways = 1;
            for (int i = 1; i < k; ++i) {
                int v = cycle[(p_idx + i) % k];
                Result res = HashWays(v, u);
                seq_hashes.push_back(res.hash_id);
                ways = (ways * res.ways) % MOD;
            }
            bool is_palin = true;
            for (int i = 0; i < (k - 1) / 2; ++i) {
                if (seq_hashes[i] != seq_hashes[k - 2 - i]) {
                    is_palin = false;
                    break;
                }
            }
            if (is_palin) ways = (ways * 2) % MOD;
            
            vector<int> rev_hashes = seq_hashes;
            reverse(rev_hashes.begin(), rev_hashes.end());
            vector<int> min_seq = min(seq_hashes, rev_hashes);
            int my_hash = get_id(2, min_seq);
            return {my_hash, ways};
        } else {
            vector<int> seq_hashes;
            long long ways = 1;
            for (int i = 0; i < k; ++i) {
                int v = cycle[i];
                Result res = HashWays(v, u);
                seq_hashes.push_back(res.hash_id);
                ways = (ways * res.ways) % MOD;
            }
            
            int sym = 0;
            vector<int> rev_hashes = seq_hashes;
            reverse(rev_hashes.begin(), rev_hashes.end());
            
            for (int shift = 0; shift < k; ++shift) {
                bool match = true;
                for (int i = 0; i < k; ++i) {
                    if (seq_hashes[(i + shift) % k] != seq_hashes[i]) {
                        match = false; break;
                    }
                }
                if (match) sym++;
                
                match = true;
                for (int i = 0; i < k; ++i) {
                    if (rev_hashes[(i + shift) % k] != seq_hashes[i]) {
                        match = false; break;
                    }
                }
                if (match) sym++;
            }
            ways = (ways * sym) % MOD;
            
            vector<int> min_seq = seq_hashes;
            for (int shift = 0; shift < k; ++shift) {
                vector<int> temp(k);
                for (int i = 0; i < k; ++i) temp[i] = seq_hashes[(i + shift) % k];
                min_seq = min(min_seq, temp);
                
                for (int i = 0; i < k; ++i) temp[i] = rev_hashes[(i + shift) % k];
                min_seq = min(min_seq, temp);
            }
            int my_hash = get_id(2, min_seq);
            return {my_hash, ways};
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N >> M)) return 0;
    
    fact[0] = 1;
    for (int i = 1; i <= 400; ++i) {
        fact[i] = (fact[i - 1] * i) % MOD;
    }
    
    adj.resize(N + 1);
    for (int i = 0; i < M; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    if (N == 1) {
        cout << 1 << "\n";
        return 0;
    }
    
    dfs(1, 0, 0);
    
    for (auto& cyc : cycles) {
        int k = cyc.size();
        for (int i = 0; i < k; ++i) {
            int u = cyc[i] + 1;
            int v = cyc[(i + 1) % k] + 1;
            cycle_edges.insert({min(u, v), max(u, v)});
        }
    }
    
    for (int i = 1; i <= N; ++i) {
        for (int v : adj[i]) {
            if (i < v) {
                if (cycle_edges.find({i, v}) == cycle_edges.end()) {
                    bridges.push_back({i, v});
                }
            }
        }
    }
    
    V_cnt = N;
    C_cnt = cycles.size();
    B_cnt = bridges.size();
    
    adj_T.resize(V_cnt + C_cnt + B_cnt);
    for (int i = 0; i < C_cnt; ++i) {
        int c_idx = N + i;
        for (int u : cycles[i]) {
            adj_T[c_idx].push_back(u);
            adj_T[u].push_back(c_idx);
        }
    }
    
    for (int i = 0; i < B_cnt; ++i) {
        int b_idx = N + C_cnt + i;
        int u = bridges[i].first - 1;
        int v = bridges[i].second - 1;
        adj_T[b_idx].push_back(u);
        adj_T[u].push_back(b_idx);
        adj_T[b_idx].push_back(v);
        adj_T[v].push_back(b_idx);
    }
    
    int A = bfs(0);
    int B = bfs(A);
    
    vector<int> path;
    int curr = B;
    while (curr != -1) {
        path.push_back(curr);
        curr = par_T[curr];
    }
    
    int center = path[path.size() / 2];
    
    Result res = HashWays(center, -1);
    cout << res.ways << "\n";
    
    return 0;
}