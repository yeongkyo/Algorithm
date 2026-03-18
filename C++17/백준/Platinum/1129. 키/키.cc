#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int N;
vector<int> H;
vector<bool> used;
vector<int> ans;
vector<int> current_path;
int D;
bool found = false;

bool check_possible(int u, int start, const vector<int>& unvisited) {
    if (unvisited.empty()) return true;

    vector<int> V = unvisited;
    V.push_back(u);
    if (u != start) V.push_back(start);

    int sz = V.size();
    vector<int> deg(sz, 0);

    for (int i = 0; i < sz; ++i) {
        for (int j = i + 1; j < sz; ++j) {
            if (abs(H[V[i]] - H[V[j]]) <= D) {
                deg[i]++;
                deg[j]++;
            }
        }
    }

    vector<bool> vis(sz, false);
    int head = 0, tail = 0;
    vector<int> q(sz);
    q[tail++] = 0;
    vis[0] = true;
    int connected_count = 1;

    while (head < tail) {
        int curr = q[head++];
        for (int i = 0; i < sz; ++i) {
            if (!vis[i] && abs(H[V[curr]] - H[V[i]]) <= D) {
                vis[i] = true;
                q[tail++] = i;
                connected_count++;
            }
        }
    }

    if (connected_count != sz) return false;

    for (int i = 0; i < sz; ++i) {
        if (V[i] == u || V[i] == start) {
            if (u == start) {
                if (deg[i] < 2 && sz > 1) return false;
            } else {
                if (deg[i] < 1) return false;
            }
        } else {
            if (deg[i] < 2) return false;
        }
    }
    return true;
}

void dfs(int u) {
    if (found) return;
    if (current_path.size() == N) {
        if (abs(H[current_path.back()] - H[current_path[0]]) <= D) {
            ans = current_path;
            found = true;
        }
        return;
    }

    vector<int> unvisited;
    for (int i = 0; i < N; ++i) {
        if (!used[i]) unvisited.push_back(i);
    }

    if (!check_possible(u, current_path[0], unvisited)) {
        return;
    }

    int last_val = -1;
    for (int i = 0; i < N; ++i) {
        if (!used[i]) {
            if (H[i] == last_val) continue;
            last_val = H[i];
            if (abs(H[i] - H[u]) <= D) {
                used[i] = true;
                current_path.push_back(i);
                dfs(i);
                if (found) return;
                current_path.pop_back();
                used[i] = false;
            }
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N)) return 0;

    H.resize(N);
    for (int i = 0; i < N; ++i) {
        cin >> H[i];
    }

    sort(H.begin(), H.end());

    D = 0;
    if (N >= 3) {
        for (int i = 0; i < N - 2; ++i) {
            D = max(D, H[i+2] - H[i]);
        }
    }

    while (!found) {
        used.assign(N, false);
        used[0] = true;
        current_path.clear();
        current_path.push_back(0);
        dfs(0);
        if (!found) {
            D++;
        }
    }

    for (int i = 0; i < N; ++i) {
        cout << H[ans[i]] << (i == N - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}