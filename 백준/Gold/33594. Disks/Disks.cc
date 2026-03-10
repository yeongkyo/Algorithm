#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef long long ll;

struct Disk {
    ll x, y, r;
};

int n;
vector<Disk> disks;
vector<vector<int>> adj;
vector<int> color;
int c1, c2;
bool isBipartite;

void dfs(int u, int c) {
    color[u] = c;
    if (c == 1) c1++;
    else c2++;

    for (int v : adj[u]) {
        if (color[v] == 0) {
            dfs(v, 3 - c);
        } else if (color[v] == c) {
            isBipartite = false;
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;
    disks.resize(n);
    adj.resize(n);
    color.assign(n, 0);

    for (int i = 0; i < n; i++) {
        cin >> disks[i].x >> disks[i].y >> disks[i].r;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            ll dx = disks[i].x - disks[j].x;
            ll dy = disks[i].y - disks[j].y;
            ll sr = disks[i].r + disks[j].r;
            if (dx * dx + dy * dy == sr * sr) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    bool possible = false;
    for (int i = 0; i < n; i++) {
        if (color[i] == 0) {
            isBipartite = true;
            c1 = 0; c2 = 0;
            dfs(i, 1);
            if (isBipartite && c1 != c2) {
                possible = true;
                break;
            }
        }
    }

    if (possible) cout << "YES" << endl;
    else cout << "NO" << endl;

    return 0;
}