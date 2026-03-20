#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

struct Region {
    long long lx, ly, hx, hy;
    int t;
};

const long long INF = 2e18;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long a, b, c, d;
    if (!(cin >> a >> b >> c >> d)) return 0;

    int n;
    cin >> n;

    vector<Region> regions(n);
    vector<long long> Vx, Vy;
    Vx.push_back(a);
    Vx.push_back(c);
    Vy.push_back(b);
    Vy.push_back(d);

    for (int i = 0; i < n; ++i) {
        cin >> regions[i].lx >> regions[i].ly >> regions[i].hx >> regions[i].hy >> regions[i].t;
        Vx.push_back(regions[i].lx);
        Vx.push_back(regions[i].hx);
        Vy.push_back(regions[i].ly);
        Vy.push_back(regions[i].hy);
    }

    sort(Vx.begin(), Vx.end());
    Vx.erase(unique(Vx.begin(), Vx.end()), Vx.end());

    sort(Vy.begin(), Vy.end());
    Vy.erase(unique(Vy.begin(), Vy.end()), Vy.end());

    int Nx = Vx.size();
    int Ny = Vy.size();
    int cell_Ny = Ny - 1;

    vector<int> W_cell;
    if (Nx > 1 && Ny > 1) {
        W_cell.assign((Nx - 1) * cell_Ny, 10);
        for (int k = 0; k < n; ++k) {
            int x1 = lower_bound(Vx.begin(), Vx.end(), regions[k].lx) - Vx.begin();
            int x2 = lower_bound(Vx.begin(), Vx.end(), regions[k].hx) - Vx.begin();
            int y1 = lower_bound(Vy.begin(), Vy.end(), regions[k].ly) - Vy.begin();
            int y2 = lower_bound(Vy.begin(), Vy.end(), regions[k].hy) - Vy.begin();

            for (int i = x1; i < x2; ++i) {
                for (int j = y1; j < y2; ++j) {
                    W_cell[i * cell_Ny + j] = regions[k].t;
                }
            }
        }
    }

    vector<long long> dist(Nx * Ny, INF);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;

    int start_x = lower_bound(Vx.begin(), Vx.end(), a) - Vx.begin();
    int start_y = lower_bound(Vy.begin(), Vy.end(), b) - Vy.begin();
    int end_x = lower_bound(Vx.begin(), Vx.end(), c) - Vx.begin();
    int end_y = lower_bound(Vy.begin(), Vy.end(), d) - Vy.begin();

    int start_node = start_x * Ny + start_y;
    dist[start_node] = 0;
    pq.push({0, start_node});

    while (!pq.empty()) {
        auto [d_u, u] = pq.top();
        pq.pop();

        int ux = u / Ny;
        int uy = u % Ny;

        if (d_u > dist[u]) continue;

        if (ux == end_x && uy == end_y) {
            cout << d_u << "\n";
            return 0;
        }

        if (ux > 0) {
            long long length = Vx[ux] - Vx[ux - 1];
            long long cost_per_unit = 10;
            if (uy > 0 && uy < Ny - 1) {
                int w1 = W_cell[(ux - 1) * cell_Ny + uy - 1];
                int w2 = W_cell[(ux - 1) * cell_Ny + uy];
                if (w1 == w2 && w1 != 10) cost_per_unit = w1;
            }
            long long new_d = d_u + length * cost_per_unit;
            int next_u = (ux - 1) * Ny + uy;
            if (new_d < dist[next_u]) {
                dist[next_u] = new_d;
                pq.push({new_d, next_u});
            }
        }

        if (ux < Nx - 1) {
            long long length = Vx[ux + 1] - Vx[ux];
            long long cost_per_unit = 10;
            if (uy > 0 && uy < Ny - 1) {
                int w1 = W_cell[ux * cell_Ny + uy - 1];
                int w2 = W_cell[ux * cell_Ny + uy];
                if (w1 == w2 && w1 != 10) cost_per_unit = w1;
            }
            long long new_d = d_u + length * cost_per_unit;
            int next_u = (ux + 1) * Ny + uy;
            if (new_d < dist[next_u]) {
                dist[next_u] = new_d;
                pq.push({new_d, next_u});
            }
        }

        if (uy > 0) {
            long long length = Vy[uy] - Vy[uy - 1];
            long long cost_per_unit = 10;
            if (ux > 0 && ux < Nx - 1) {
                int w1 = W_cell[(ux - 1) * cell_Ny + uy - 1];
                int w2 = W_cell[ux * cell_Ny + uy - 1];
                if (w1 == w2 && w1 != 10) cost_per_unit = w1;
            }
            long long new_d = d_u + length * cost_per_unit;
            int next_u = ux * Ny + uy - 1;
            if (new_d < dist[next_u]) {
                dist[next_u] = new_d;
                pq.push({new_d, next_u});
            }
        }

        if (uy < Ny - 1) {
            long long length = Vy[uy + 1] - Vy[uy];
            long long cost_per_unit = 10;
            if (ux > 0 && ux < Nx - 1) {
                int w1 = W_cell[(ux - 1) * cell_Ny + uy];
                int w2 = W_cell[ux * cell_Ny + uy];
                if (w1 == w2 && w1 != 10) cost_per_unit = w1;
            }
            long long new_d = d_u + length * cost_per_unit;
            int next_u = ux * Ny + uy + 1;
            if (new_d < dist[next_u]) {
                dist[next_u] = new_d;
                pq.push({new_d, next_u});
            }
        }
    }

    return 0;
}