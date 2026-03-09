#include <iostream>
#include <vector>

using namespace std;

int h, w;
char grid[10][10];
bool is_X[10][10];
vector<pair<int, int>> cells;
int N;

bool found_solution = false;
vector<pair<int, int>> current_deltas;
bool used_cell[10][10];

void solve_k(int k, int S, int idx_start, int num_picked) {
    if (found_solution) return;
    if (num_picked == k - 1) {
        pair<int, int> p1_cells[100];
        int p1_size = 0;
        bool valid = true;
        
        pair<int, int> modified[500];
        int mod_cnt = 0;
        
        for (int j = 0; j < N; ++j) {
            int r = cells[j].first;
            int c = cells[j].second;
            if (used_cell[r][c]) continue;
            
            if (p1_size == S) { valid = false; break; }
            p1_cells[p1_size++] = {r, c};
            
            for (int i = 0; i < k; ++i) {
                int nr = r + current_deltas[i].first;
                int nc = c + current_deltas[i].second;
                if (nr < 0 || nr >= h || nc < 0 || nc >= w || !is_X[nr][nc] || used_cell[nr][nc]) {
                    valid = false;
                    break;
                }
                used_cell[nr][nc] = true;
                modified[mod_cnt++] = {nr, nc};
            }
            if (!valid) break;
        }
        
        for (int m = 0; m < mod_cnt; ++m) {
            used_cell[modified[m].first][modified[m].second] = false;
        }
        
        if (valid && p1_size == S) {
            int connected_count = 1;
            bool visited_p1[100] = {false};
            visited_p1[0] = true;
            int q[100];
            int head = 0, tail = 0;
            q[tail++] = 0;
            
            while (head < tail) {
                int u = q[head++];
                for (int v = 0; v < S; ++v) {
                    if (!visited_p1[v]) {
                        int dr = p1_cells[u].first - p1_cells[v].first;
                        int dc = p1_cells[u].second - p1_cells[v].second;
                        if (dr * dr + dc * dc == 1) {
                            visited_p1[v] = true;
                            q[tail++] = v;
                            connected_count++;
                        }
                    }
                }
            }
            
            if (connected_count == S) {
                found_solution = true;
                char ans[10][10];
                for (int r = 0; r < h; ++r) {
                    for (int c = 0; c < w; ++c) {
                        ans[r][c] = '.';
                    }
                }
                
                for (int j = 0; j < N; ++j) {
                    int r = cells[j].first;
                    int c = cells[j].second;
                    if (used_cell[r][c]) continue;
                    
                    for (int i = 0; i < k; ++i) {
                        int nr = r + current_deltas[i].first;
                        int nc = c + current_deltas[i].second;
                        ans[nr][nc] = '1' + i;
                        used_cell[nr][nc] = true;
                    }
                }
                
                for (int r = 0; r < h; ++r) {
                    for (int c = 0; c < w; ++c) {
                        cout << ans[r][c];
                    }
                    cout << "\n";
                }
            }
        }
        return;
    }
    
    for (int i = idx_start; i < N; ++i) {
        int dr = cells[i].first - cells[0].first;
        int dc = cells[i].second - cells[0].second;
        current_deltas.push_back({dr, dc});
        solve_k(k, S, i + 1, num_picked + 1);
        current_deltas.pop_back();
        if (found_solution) return;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> h >> w)) return 0;
    
    for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
            cin >> grid[r][c];
            if (grid[r][c] == 'X') {
                is_X[r][c] = true;
                cells.push_back({r, c});
            } else {
                is_X[r][c] = false;
            }
        }
    }
    
    N = cells.size();
    
    for (int k = 2; k <= 5; ++k) {
        if (N % k != 0) continue;
        int S = N / k;
        
        current_deltas.clear();
        current_deltas.push_back({0, 0});
        solve_k(k, S, 1, 0);
        
        if (found_solution) break;
    }
    
    if (!found_solution) {
        cout << "No solution\n";
    }
    
    return 0;
}