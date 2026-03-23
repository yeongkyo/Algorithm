#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

struct TilePos {
    int r, c1, c2;
};

int n;
int tileID[505][1005];
int tileVal[505][1005];
TilePos posMap[250005];
int parent[250005];
int dist[250005];

int dr[] = {-1, 1, 0, 0};
int dc[] = {0, 0, -1, 1};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    int total_tiles = n * n - n / 2;
    int currentID = 1;
    for (int r = 1; r <= n; r++) {
        if (r % 2 != 0) {
            for (int i = 1; i <= n; i++) {
                int a, b;
                cin >> a >> b;
                int c1 = 2 * i - 1;
                int c2 = 2 * i;
                tileID[r][c1] = tileID[r][c2] = currentID;
                tileVal[r][c1] = a;
                tileVal[r][c2] = b;
                posMap[currentID] = {r, c1, c2};
                currentID++;
            }
        } else {
            for (int i = 1; i <= n - 1; i++) {
                int a, b;
                cin >> a >> b;
                int c1 = 2 * i;
                int c2 = 2 * i + 1;
                tileID[r][c1] = tileID[r][c2] = currentID;
                tileVal[r][c1] = a;
                tileVal[r][c2] = b;
                posMap[currentID] = {r, c1, c2};
                currentID++;
            }
        }
    }

    for (int i = 1; i <= total_tiles; i++) dist[i] = -1;

    queue<int> q;
    q.push(1);
    dist[1] = 1;
    int maxID = 1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u > maxID) maxID = u;

        int r = posMap[u].r;
        int cols[2] = {posMap[u].c1, posMap[u].c2};

        for (int i = 0; i < 2; i++) {
            int c = cols[i];
            for (int d = 0; d < 4; d++) {
                int nr = r + dr[d];
                int nc = c + dc[d];

                if (nr < 1 || nr > n || nc < 1 || nc > 2 * n) continue;
                int v = tileID[nr][nc];
                if (v != 0 && v != u && dist[v] == -1) {
                    if (tileVal[r][c] == tileVal[nr][nc]) {
                        dist[v] = dist[u] + 1;
                        parent[v] = u;
                        q.push(v);
                    }
                }
            }
        }
    }

    vector<int> path;
    int curr = maxID;
    while (curr != 0) {
        path.push_back(curr);
        curr = parent[curr];
    }
    reverse(path.begin(), path.end());

    cout << path.size() << "\n";
    for (int i = 0; i < path.size(); i++) {
        cout << path[i] << (i == path.size() - 1 ? "" : " ");
    }
    cout << "\n";

    return 0;
}