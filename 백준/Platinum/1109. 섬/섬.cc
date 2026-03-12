#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

int N, M;
string grid[55];
int island_id[55][55];
int sea_id[55][55];
int num_islands = 0;
int num_seas = 0;

int dr4[] = {-1, 1, 0, 0};
int dc4[] = {0, 0, -1, 1};
int dr8[] = {-1, -1, -1, 0, 0, 1, 1, 1};
int dc8[] = {-1, 0, 1, -1, 1, -1, 0, 1};

vector<int> island_children[3005]; 
vector<int> sea_children[3005];    
int height[3005];

int get_island_height(int u) {
    if (height[u] != -1) return height[u];
    int max_child_h = -1;
    for (int sea : island_children[u]) {
        for (int child_island : sea_children[sea]) {
            max_child_h = max(max_child_h, get_island_height(child_island));
        }
    }
    height[u] = max_child_h + 1;
    return height[u];
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N >> M)) return 0;

    // 격자를 바다('.')로 한 겹 패딩
    grid[0] = string(M + 2, '.');
    for (int i = 1; i <= N; ++i) {
        string s;
        cin >> s;
        grid[i] = "." + s + ".";
    }
    grid[N + 1] = string(M + 2, '.');

    // 1. 섬 연결 요소 파악 (8방향)
    for (int r = 0; r < N + 2; ++r) {
        for (int c = 0; c < M + 2; ++c) {
            if (grid[r][c] == 'x' && island_id[r][c] == 0) {
                num_islands++;
                queue<pair<int, int>> q;
                q.push({r, c});
                island_id[r][c] = num_islands;
                while (!q.empty()) {
                    auto [cr, cc] = q.front(); q.pop();
                    for (int i = 0; i < 8; ++i) {
                        int nr = cr + dr8[i];
                        int nc = cc + dc8[i];
                        if (nr >= 0 && nr < N + 2 && nc >= 0 && nc < M + 2) {
                            if (grid[nr][nc] == 'x' && island_id[nr][nc] == 0) {
                                island_id[nr][nc] = num_islands;
                                q.push({nr, nc});
                            }
                        }
                    }
                }
            }
        }
    }

    if (num_islands == 0) {
        cout << -1 << "\n";
        return 0;
    }

    // 2. 바다 연결 요소 파악 (4방향)
    for (int r = 0; r < N + 2; ++r) {
        for (int c = 0; c < M + 2; ++c) {
            if (grid[r][c] == '.' && sea_id[r][c] == 0) {
                num_seas++;
                queue<pair<int, int>> q;
                q.push({r, c});
                sea_id[r][c] = num_seas;
                while (!q.empty()) {
                    auto [cr, cc] = q.front(); q.pop();
                    for (int i = 0; i < 4; ++i) {
                        int nr = cr + dr4[i];
                        int nc = cc + dc4[i];
                        if (nr >= 0 && nr < N + 2 && nc >= 0 && nc < M + 2) {
                            if (grid[nr][nc] == '.' && sea_id[nr][nc] == 0) {
                                sea_id[nr][nc] = num_seas;
                                q.push({nr, nc});
                            }
                        }
                    }
                }
            }
        }
    }

    // 각 그룹별 좌표 모아두기
    vector<vector<pair<int,int>>> island_cells(num_islands + 1);
    vector<vector<pair<int,int>>> sea_cells(num_seas + 1);
    for (int r = 0; r < N + 2; ++r) {
        for (int c = 0; c < M + 2; ++c) {
            if (island_id[r][c] > 0) island_cells[island_id[r][c]].push_back({r, c});
            if (sea_id[r][c] > 0) sea_cells[sea_id[r][c]].push_back({r, c});
        }
    }

    vector<bool> vis_island(num_islands + 1, false);
    vector<bool> vis_sea(num_seas + 1, false);

    // 가장 바깥쪽 바다를 루트로 트리 구성
    int out_sea = sea_id[0][0];
    queue<pair<int, int>> tree_q;
    tree_q.push({0, out_sea});
    vis_sea[out_sea] = true;

    while (!tree_q.empty()) {
        auto [type, id] = tree_q.front();
        tree_q.pop();

        if (type == 0) {
            for (auto [r, c] : sea_cells[id]) {
                for (int i = 0; i < 8; ++i) {
                    int nr = r + dr8[i];
                    int nc = c + dc8[i];
                    if (nr >= 0 && nr < N + 2 && nc >= 0 && nc < M + 2) {
                        if (grid[nr][nc] == 'x') {
                            int nxt = island_id[nr][nc];
                            if (!vis_island[nxt]) {
                                vis_island[nxt] = true;
                                sea_children[id].push_back(nxt);
                                tree_q.push({1, nxt});
                            }
                        }
                    }
                }
            }
        } else {
            for (auto [r, c] : island_cells[id]) {
                for (int i = 0; i < 8; ++i) {
                    int nr = r + dr8[i];
                    int nc = c + dc8[i];
                    if (nr >= 0 && nr < N + 2 && nc >= 0 && nc < M + 2) {
                        if (grid[nr][nc] == '.') {
                            int nxt = sea_id[nr][nc];
                            if (!vis_sea[nxt]) {
                                vis_sea[nxt] = true;
                                island_children[id].push_back(nxt);
                                tree_q.push({0, nxt});
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 1; i <= num_islands; ++i) height[i] = -1;

    int max_h = -1;
    for (int i = 1; i <= num_islands; ++i) {
        if (height[i] == -1) {
            get_island_height(i);
        }
        max_h = max(max_h, height[i]);
    }

    vector<int> cnt(max_h + 1, 0);
    for (int i = 1; i <= num_islands; ++i) {
        cnt[height[i]]++;
    }

    for (int i = 0; i <= max_h; ++i) {
        cout << cnt[i] << (i == max_h ? "" : " ");
    }
    cout << "\n";

    return 0;
}